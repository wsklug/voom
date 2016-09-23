// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   HoHai Van and William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
// 
// 
//----------------------------------------------------------------------

#include "NobleElem.h"
#include "math.h"
#include "ShapeQ4.h"
#include "QuadQuadrature.h"
#include <vector>
#include <iostream>

namespace voom
{

  NobleElem::QuadPointStruct::
  QuadPointStruct(double w, const NobleElem::Shape_t & s, NodeMatrix_t & nds) {
    weight=w;
    shapeFunctions=s.functions();

    // Compute spatial derivatives of shape functions from
    // parametric derivatives by transforming with matrix jacobian
    // of isoparametric mapping.
    
    // parametric derivatives from shape object
    const Shape_t::DerivativeContainer & dnds = s.derivatives();

    // for each quad point:  compute the jacobian matrix 
      tvmet::Matrix<double,2,2> Jac (0.0);
            
    // for each quadrature point:  loop over the nodes of the element 
      for (int i = 0; i <= 3 ; i++){
        
        double xi=nds(i,0);
        double yi=nds(i,1);
        
 
	Jac(0,0)+=dnds[i](0)*xi;
	Jac(0,1)+=dnds[i](1)*xi;
	Jac(1,0)+=dnds[i](0)*yi;
	Jac(1,1)+=dnds[i](1)*yi;
//        cout << "X,Y  "  << xi <<","<<yi << "  dnds[i]" << dnds[i](0)<<","<<dnds[i][1]<<endl;

//      cout << "Jac  " << Jac(0,0) << "  " << Jac(0,1) << " " << Jac (1,0) << "  " << Jac (1,1);
//      cout << endl;+



      }

//      cout << "Jac  " << Jac(0,0) << "  " << Jac(0,1) << " " << Jac (1,0) << "  " << Jac (1,1);
//      cout << endl;

    // for each quadrature point:  calculate the determinant of the jacobian

      double detJ;
      detJ=Jac(0,0)*Jac(1,1)-Jac(0,1)*Jac(1,0);
      
//      cout << "detJ  "  << detJ << endl;
   
    // multiply all the weights by the determinate of the jacobian so only have to pass the weight
      weight *= detJ;

    // for each quadrature point:  calculate the inverse of the jacobian.

      tvmet::Matrix<double,2,2> IJac (0.0);

      IJac(0,0)=Jac(1,1);
      IJac(0,1)=-Jac(0,1);
      IJac(1,0)=-Jac(1,0);
      IJac(1,1)=Jac(0,0);

      IJac/=detJ;

//      cout << "IJac   " << IJac(0,0) << "  " << IJac(0,1) << " " << IJac (1,0) << "  "<<
//         IJac (1,1) << endl;


    // calculate shapeDerivatives which are dndx
      shapeDerivatives.resize(4);
      for (int i=0;i<=3;i++) {
      
        shapeDerivatives [i][0]=0.0;
        shapeDerivatives [i][1]=0.0;
      }

      for (int i=0;i<=3;i++) {
        shapeDerivatives[i](0)=IJac(0,0)*dnds[i](0)+IJac(1,0)*dnds[i](1);
        shapeDerivatives[i](1)=IJac(0,1)*dnds[i](0)+IJac(1,1)*dnds[i](1);
      }

   return;
  }


  NobleElem::NobleElem( const Quadrature_t & quad,
		        NodeMatrix_t & nodes,
                        QuadMatrix_t & quads, const double i_stim)
  {
       
    //! initialize quad points
    _quadPoints.clear();
    for(Quadrature_t::ConstPointIterator p=quad.begin(); p!=quad.end(); p++) {
      Shape_t shp( p-> coords );
      _quadPoints.push_back( QuadPointStruct(p->weight, shp, nodes) ); 

    }
    compute (i_stim, nodes, quads);
    return;
  }

  void NobleElem::compute(double i_stim, NodeMatrix_t & nodes, QuadMatrix_t & quads)
  {

    int counter=0;
    tvmet::Vector<double,4> i_ion_vector(0.0);
    tvmet::Matrix<double,4,4> s(0.0);
    tvmet::Vector<double,4> capacitance(0.0);

    // loop for every quadrature point
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){

    // get the shape functions and the derivatives
     
        const Shape_t::DerivativeContainer  & DN = p->shapeDerivatives;
        const Shape_t::FunctionContainer & N = p->shapeFunctions;

   // Quadrature of element stiffness and ionic current 
   // Already have m,n,h at the quadrature points
        double i_ion,i_na,i_k,i_leak;
        double volt_quad,m,n,h;
      
        volt_quad=quads(counter,0);
        m=quads(counter,1);
        h=quads(counter,2);
        n=quads(counter,3);

      
        i_na=(400.0*pow(m,3.0)*h+0.14)*(volt_quad-40.0);

        double g_k1,g_k2;
        g_k1=1.2*exp((-volt_quad-90.0)/50.0)+0.015*exp((volt_quad+90.0)/60.0);
        g_k2=1.2*pow(n,4.0);

        i_k=(g_k1+g_k2)*(volt_quad+100.0);

        i_leak=0.0*(volt_quad-(-60.0));

        i_ion=-(i_na+i_k+i_leak-i_stim);

        double Dx=0.00025;
        double Dy=0.00025;


        for (int i=0;i<=3;i++) {

//Calculate the ionic current vector         
          i_ion_vector(i)+=(p->weight)*N[i]*i_ion;

//Caclulate the lumped capacitance vector
          double c_m=12.0;
          capacitance(i)+=(p->weight)*N[i]*c_m;

          for (int j=0;j<=3;j++) {
//Calculate the stiffness matrix	   
      	      s(i,j)=s(i,j)+(p->weight)*(DN[i](0)*Dx*DN[j](0)
                  +DN[i](1)*Dy*DN[j](1));
             } 
       

           }
      

        counter+=1;
  }// end loop for all the quadrature points;


// calculate the dvdt vector for the element 
      tvmet::Vector<double,4> dvdt (0.0);

      for (int i=0;i<=3;i++) {
        nodes(i,4)=capacitance(i);
        dvdt(i)=i_ion_vector(i);

        for (int j=0;j<=3;j++) {
           dvdt(i)-=s(i,j)*nodes(j,2);
           } 
        nodes(i,3)=dvdt(i);
      }




}// end function compute of NobleElem

} // namespace voom
