// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    HoHai Van and William S. Klug
//                University of California Los Angeles
//                    (C) 2006 All Rights Reserved
//
// Revision 2:  8/17/06:  
// This version simulates the Beeler-Reuter model.  Changes were also made
// in the compute function to split up the voltage calculation and the time
// stepping of the gating variables.  
//
// Revision 3:  8/26/06:
// The gating variables formulations were changed to the analytical method
// which allowed the time step to be larger.  Also, the code now includes
// adaptive time stepping to update the gating variables only when dv/dt
// is a certain number.
//
// Revision 4:  9/3/06: 
// Moved the calculation of the gating variables to a lookup table called
// gate to speed up the program.  Update: did not speed up program
//
// Revision 5: 9 9/5/06:
// Changed to Luo-Rudy model
//----------------------------------------------------------------------

#include "CardiacPotential.h"
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <blitz/array.h>
#include "math.h"

namespace voom {

  /*!  The QuadPointStruct constructor has to initialize the spatial
    derivatives of the shape functions as well as the effective
    quadrature weight (scaled by the scalar jacobian).  First we     
    compute the Jacobian
    \f[
    J = \det(J_{i\alpha}), \qquad J_{i\alpha} = \sum_{a=1}^{\cal N} N_{a,\alpha} x_{ia} , \qquad N_{a,\alpha} = \frac{\partial N_a}{\partial s^\alpha}
    \f]
    and then compute the spatial derivatives of the shape functions
    \f[
    \frac{\partial N_a}{\partial x_i} = J^{-1}_{\alpha i} \frac{\partial N_a}{\partial s_\alpha}
    \f]
   */
//  template< class Quadrature_t,
//	    class Shape_t >
  CardiacPotential::QuadPointStruct::
  QuadPointStruct(double w, 
		  const Shape_t & s, 
		  const CardiacPotential::VoltageNodeContainer & nds)
  {
    weight = w;
    shapeFunctions = s.functions();
    
    // Compute spatial derivatives of shape functions from
    // parametric derivatives by transforming with matrix jacobian
    // of isoparametric mapping.
    
    // parametric derivatives from shape object
    const Shape_t::DerivativeContainer & dnds = s.derivatives();

    // matrix jacobian dxds;
    tvmet::Matrix<double,2,2> dxds(0.0);
    for(int i=0; i<2; i++) {
      for(int alpha=0; alpha<2; alpha++) {	   
	for(int a=0; a<nds.size(); a++) {
	  dxds(i,alpha) += dnds[a](alpha)*( nds[a]->position(i) );
	}
      }	      
    }

    // compute scalar jacobian and scale the quadrature weight with it
    double J=dxds(0,0)*dxds(1,1) - dxds(0,1)*dxds(1,0);
    weight *= J;
	
    // invert matrix jacobian
    tvmet::Matrix<double,2,2> invJac(0.0);
    invJac(0,0) = dxds(1,1);
    invJac(1,1) = dxds(0,0);
    invJac(0,1) = -dxds(0,1);
    invJac(1,0) = -dxds(1,0);
    invJac /= J;

    // spatial derivatives dndx
    shapeDerivatives.resize(dnds.size());
    for(int a=0; a<nds.size(); a++) {	  
      for(int i=0; i<2; i++) {
	shapeDerivatives[a](i) = 0.0;
	for(int alpha=0; alpha<2; alpha++) {	   
	  shapeDerivatives[a](i) += dnds[a](alpha)*invJac(alpha,i);
	}
      }
    }

    internalNode = new InternalNode(0,NodeBase::DofIndexMap(7),InternalNode::Point(0.0));
    
    // Note: this will call the default constructor for the internal
    // nodes.  Thus they will need to be initialized later (setting
    // id and index members), when we want to use them with a
    // solver.

    // Each internal node will represent an ionic gating variable.  
    // For the Noble model, the gating variables are m,n,h
    return;
  }


  // Constructor
//    template< class Quadrature_t,
//	    class Shape_t >
    CardiacPotential::
    CardiacPotential( const Quadrature_t & quad,
				      const VoltageNodeContainer & nodes )
    {

    //! initialize NodeContainer
    unsigned nNodes = nodes.size();
    _vNodes = nodes;
  

//    for(ConstVoltageNodeIterator n=_vNodes.begin(); n!=_vNodes.end(); n++) 
//      _baseNodes.push_back(*n);
        
    //! initialize quad points
    _quadPoints.clear();
    for(Quadrature_t::ConstPointIterator p=quad.begin(); p!=quad.end(); p++) {
      Shape_t shp( p->coords );
      _quadPoints.push_back( QuadPointStruct(p->weight, shp, _vNodes) ); 
	if( shp.functions().size() != nNodes ) {
	  std::cout << "Number of nodes: " << nNodes << std::endl
		    << "Number of functions: " << shp.functions().size()
		    << std::endl
		    << "These should be equal." << std::endl;
	  exit(0);
	}
      
    }

    for(ConstQuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++) {
//       _baseNodes.push_back(p->internalNode);
      _iNodes.push_back(p->internalNode);
       }

    //! initialize the gating variables by default
    _stiffness.resize( _vNodes.size(), _vNodes.size() );
    _stiffness = 0.0;
    _volt_quad_old=0.0;
    _dt=0.001;

//<<<<<<< .mine
//    _stiffness = Matrix44(0.0);

    // Initialize the damping and force components
    for(int a=0; a<_vNodes.size(); a++) {
       _vNodes[a]->setDamping(0,0.0);
       _vNodes[a]->setForce(0,0.0);
       }

    for(QuadPointIterator p=_quadPoints.begin();p!=_quadPoints.end(); p++){
      const Shape_t::FunctionContainer &  N = p->shapeFunctions;
      const Shape_t::DerivativeContainer &  DN = p->shapeDerivatives;
      InternalNode * in = p->internalNode;

      //Compute voltage, capacitance, and stiffness matrix
      double volt_quad=0.0;
      double c_m=1.0;
      double Dx=0.001;
      double Dy=0.001;

      for(int a=0; a<_vNodes.size(); a++) {
	 volt_quad += N[a]*(_vNodes[a]->getPoint(0));
         double c = N[a]*c_m *(p->weight);
         _vNodes[a]->addDamping(0,c);
         for(int b=0;b<_vNodes.size();b++) {
            _stiffness(a,b)+=(p->weight)*(DN[a](0)*Dx*DN[b](0)
                     +DN[a](1)*Dy*DN[b](1));
            }

         }

      double a_m=0.32*(volt_quad+47.13)/(1.0-exp(-0.1*(volt_quad+47.13)));
      double b_m=0.08*exp(-volt_quad/11.0);
      double m=a_m/(a_m+b_m);
      in->setPoint(0,m);

      double a_h=0.135*exp((80.0+volt_quad)/-6.8);
      double b_h=3.56*exp(0.079*volt_quad)+310000.0*exp(0.35*volt_quad);
      double h=a_h/(a_h+b_h);
      in->setPoint(1,h);

      double a_j=(-127140.0*exp(0.2444*volt_quad)-0.00003474*exp(-0.04391*volt_quad))*(volt_quad+37.78)/(1+exp(0.311*(volt_quad+79.23)));
      double b_j=0.1212*exp(-0.01052*volt_quad)/(1+exp(-0.1378*(volt_quad+40.14)));
      double j=a_j/(a_j+b_j);
      in->setPoint(2,j);

      double a_d=(0.095*exp(-0.01*(volt_quad-5.0)))/(exp(-0.072*(volt_quad-5.0))+1.0);
      double b_d=(0.07*exp(-0.017*(volt_quad+44.0)))/(exp(0.05*(volt_quad+44.0))+1.0);
      double d=a_d/(a_d+b_d);
      in->setPoint(3,d);

      double a_f=(0.012*exp(-0.008*(volt_quad+28.0)))/(exp(0.15*(volt_quad+28.0))+1.0);
      double b_f=(0.0065*exp(-0.02*(volt_quad+30.0)))/(exp(-0.2*(volt_quad+30.0))+1.0);
      double f=a_f/(a_f+b_f);
      in->setPoint(4,f);

      double a_x=0.0005*exp(0.083*(volt_quad+50.0))/(1.0+exp(0.057*(volt_quad+50.0)));
      double b_x=0.0013*exp(-0.06*(volt_quad+20.0))/(1.0+exp(-0.04*(volt_quad+20.0)));
      double x=a_x/(a_x+b_x);
      in->setPoint(5,x);

      in->setPoint(6,0.0002);  // Initial calcium concentration
      } //end loop of quad points

  return;
  }

  /*!  The compute method needs to compute the nodal ionic and
    diffusive currents as well as the (diagonal) lumped conduction
    matrix coefficients.  
    The currents are written 
    \f[
    f^{ion}_{a} = \int_\Omega N_a \frac{I_{ion}}{C} dV
    \f]
    and 
    \f[
    f^{dif}_{a} = \int_\Omega \nabla N_a \cdot (D\nabla V) dV 
    = \int_\Omega  N_{a,i} \cdot (D_{ij} V_{,j}) dV ,
    \f]
    and the row-summed conduction coefficients are 
    \f[
    C_a = \sum_{b=1}^{\cal N} C_{ab} = \sum_{b=1}^{\cal N} \int_\Omega N_a N_b dV .
    \f]
    Assuming the shape functions satisfy theaddDamping partition of unity property, 
    the conduction coefficients simplify to 
    \f[
    C_a = \int_\Omega N_a dV .
    \f]
    
    Approximating these nodal integrals by quadrature, we obtain
    \f[
    f^{ion}_{a} = \sum_{p=1}^Q (N_a \frac{I_{ion}}{C} J )|_{s^p} w_p
    \f]
    \f[
    f^{dif}_{a} = \sum_{p=1}^Q ( N_{a,i} \cdot (D_{ij} V_{,j}) J )|_{s^p} w_p  
    \f]
    \f[
    C_a = \sum_{p=1}^Q ( N_a J )|_{s^p} w_p 
    \f]
    where \f$s^p\f$ are the parametric coordinates of the Gauss points 
    and \f$J\f$ is the jacobian of the isoparametric mapping.

    We loop over the quad points, and at each quad point perform the 
    following steps:
    <ol>
    <li> Compute the voltage and its gradient
    \f[ 
    V = \sum_{a} N_a V_a \qquad V_i = \sum_a N_{a,i} V_a
    \f]
    <li> Compute the ionic current \f$I^{ion}(V)\f$
    <li> Compute quad point contributions to the sums for \f$f_a^{ion}\f$, 
    \f$f_a^{dif}\f$, and \f$C_a\f$.  Add these to the $a^{th}$ node using
    the methods updateForce, and addDamping.
    </ol>      
  */
  // for now we won't worry about the boolean f1, etc.  Just compute
  // forces and lumped conduction matrix
//<<<<<<< .mine
  void CardiacPotential::compute_v( double i_stim) 
//=======
//  template< class Quadrature_t,
//	    class Shape_t >
//  void CardiacPotential<Quadrature_t,Shape_t>::
//  compute( bool f0, double i_stim, bool f1, double dt, bool f2) 
//>>>>>>> .r284
  {
    //initialize f_ion or Force
//    for(int a=0; a<_vNodes.size(); a++) {
//          _vNodes[a]->setForce(0.0);
 //      }

    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
    
        const Shape_t::FunctionContainer &  N = p->shapeFunctions;
        const Shape_t::DerivativeContainer &  DN = p->shapeDerivatives;
        InternalNode * in = p->internalNode;
      
        //Compute voltage
        double volt_quad=0.0;
        for(int a=0; a<_vNodes.size(); a++) {
           volt_quad += N[a]*(_vNodes[a]->getPoint(0));
           }
      
        // Compute I^{ion}
        double i_ion = 0.0; 

        // Quadrature of ionic current 
        // Already have m,n,j,d,f,xi,calcium at the quadrature points  
        double m,h,j,d,f,x,ca;
      
        m= in->getPoint(0);
        h= in->getPoint(1);
        j= in->getPoint(2);
        d= in->getPoint(3);
        f= in->getPoint(4);
        x= in->getPoint(5);
        ca= in->getPoint(6);

        double i_na=23.0*pow(m,3.0)*h*j*(volt_quad-54.4);

        double xi; 
        if (volt_quad<=-100.0) {
           xi=1.0;
           }
        else {
           xi=2.837*(exp(0.04*(volt_quad+77.0))-1.0)/((volt_quad+77.0)*exp(0.04*(volt_quad+35.0)));
           }
        double gk_bar=0.282;
        double i_k=gk_bar*xi*x*(volt_quad+77.0);

        double E_si=7.7-13.0287*log(ca);
        double i_si=0.09*d*f*(volt_quad-E_si);

        double gk1_bar=0.6047;
        double a_k1=1.02/(1.0+exp(0.2385*(volt_quad+87.95-59.215)));
        double b_k1=0.49124*exp(0.08032*(volt_quad+87.95+5.476))+exp(0.06175*(volt_quad+87.95-594.31));
        b_k1/=1.0+exp(-0.5143*(volt_quad+87.95+4.753));
        double i_k1=gk1_bar*a_k1/(a_k1+b_k1)*(volt_quad+87.95);

        double kp=1.0/(1.0+exp((7.488-volt_quad)/5.98));
        double i_kp=0.0183*kp*(volt_quad+87.95);

        double i_b=0.03921*(volt_quad+59.87);

        i_ion=-(i_na+i_k+i_si+i_k1+i_kp+i_b-i_stim);

        // Compute f_a^{ion}
        double f_ion;
        for(int a=0; a<_vNodes.size(); a++) {
    	   f_ion = N[a]*i_ion*(p->weight);
           _vNodes[a]->addForce(0,f_ion);
          
           }
  
       } // end over quad points

//      for(int a=0;a<_vNodes.size();a++) {
//         double dummy2=_vNodes[a]->getForce(0);
//         std::cout<<a<<" before stiffness "<<dummy2<<std::endl;
//      }
      
      for(int a=0; a<_vNodes.size(); a++) {
          for(int b=0;b<_vNodes.size();b++) {
              double stiff=-_stiffness(a,b)*(_vNodes[b]->getPoint(0));
              _vNodes[a]->addForce(0,stiff);
             }
          
         }

//    for(int a=0;a<_vNodes.size();a++) {
//         double dummy2=_vNodes[a]->getForce(0);
//         std::cout<<a<<" after stiffness "<<dummy2<<std::endl;
//      }
     return;
     }// end compute_v

  void CardiacPotential::compute_gate() {

    // this portion of compute advances the internal nodes (gating variables) by dt
    for(QuadPointIterator p=_quadPoints.begin(); 
       p!=_quadPoints.end(); p++){
       const Shape_t::FunctionContainer &  N = p->shapeFunctions;
       InternalNode * in = p->internalNode;
       //Compute voltage
       double volt_quad=0.0;
       for(int a=0; a<_vNodes.size(); a++) {
          volt_quad += N[a]*(_vNodes[a]->getPoint(0));
          }

       if ((fabs(volt_quad-_volt_quad_old)>=0.02)||(_dt>0.0099)) {

          //reset _volt_quad_old           
          _volt_quad_old=volt_quad;

          //Calculate gating variables
          double m,h,j,d,f,x,ca;
      
          m= in->getPoint(0);
          h= in->getPoint(1);
          j= in->getPoint(2);
          d= in->getPoint(3);
          f= in->getPoint(4);
          x= in->getPoint(5);
          ca= in->getPoint(6);
          
          double dt=_dt;

          double E_si=7.7-13.0287*log(ca);
          double i_si=0.09*d*f*(volt_quad-E_si);
          ca+=dt*(-0.0001*i_si+0.07*(0.0001-ca));
          in->setPoint(6,ca);
          
          double a_m=0.32*(volt_quad+47.13)/(1.0-exp(-0.1*(volt_quad+47.13)));
          double b_m=0.08*exp(-volt_quad/11.0);
          double m_tau=1.0/(a_m+b_m);
          double m_inf=a_m*m_tau;
          m=m_inf-(m_inf-m)*exp(-dt/m_tau);
          in->setPoint(0,m);
          
          double a_h,b_h,a_j,b_j;
          if (volt_quad>=-40.0) {
             a_h=0.0;
             b_h=1.0/(0.13*(1.0+exp((volt_quad+10.66)/-11.1)));
             a_j=0.0;
             b_j=0.3*exp(-0.0000002535*volt_quad)/(1.0+exp(-0.1*(volt_quad+32.0)));
             }
          else {
             a_h=0.135*exp((80.0+volt_quad)/-6.8);
             b_h=3.56*exp(0.079*volt_quad)+310000.0*exp(0.35*volt_quad);
             a_j=(-127140.0*exp(0.2444*volt_quad)-0.00003474*exp(-0.04391*volt_quad))*(volt_quad+37.78)/(1+exp(0.311*(volt_quad+79.23)));
             b_j=0.1212*exp(-0.01052*volt_quad)/(1.0+exp(-0.1378*(volt_quad+40.14)));

             }  

          double h_tau=1.0/(a_h+b_h);
          double h_inf=a_h*h_tau;
          h=h_inf-(h_inf-h)*exp(-dt/h_tau);
          in->setPoint(1,h);

          double j_tau=1.0/(a_j+b_j);
          double j_inf=a_j*j_tau;
          j=j_inf-(j_inf-j)*exp(-dt/j_tau);
          in->setPoint(2,j);

          double a_d=(0.095*exp(-0.01*(volt_quad-5.0)))/(exp(-0.072*(volt_quad-5.0))+1.0);
          double b_d=(0.07*exp(-0.017*(volt_quad+44.0)))/(exp(0.05*(volt_quad+44.0))+1.0);
          double d_tau=1.0/(a_d+b_d);
          double d_inf=a_d*d_tau;
          d=d_inf-(d_inf-d)*exp(-dt/d_tau);
          in->setPoint(3,d);

          double a_f=(0.012*exp(-0.008*(volt_quad+28.0)))/(exp(0.15*(volt_quad+28.0))+1.0);
          double b_f=(0.0065*exp(-0.02*(volt_quad+30.0)))/(exp(-0.2*(volt_quad+30.0))+1.0);
          double f_tau=1.0/(a_f+b_f);
          double f_inf=a_f*f_tau;
          f=f_inf-(f_inf-f)*exp(-dt/f_tau);
          in->setPoint(4,f);

          double a_x=0.0005*exp(0.083*(volt_quad+50.0))/(1.0+exp(0.057*(volt_quad+50.0)));
          double b_x=0.0013*exp(-0.06*(volt_quad+20.0))/(1.0+exp(-0.04*(volt_quad+20.0)));
          double x_tau=1.0/(a_x+b_x);
          double x_inf=a_x*x_tau;
          x=x_inf-(x_inf-x)*exp(-dt/x_tau);
          in->setPoint(5,x);

          _dt=0.001; //reset _dt
          }//end if statement
       // else small change in dv so do not update gating variables

       else {
          _dt+=0.001;
          }


      } //end loop of quad points
    return;
    }  //end compute_gate

  
}; // namespace voom
