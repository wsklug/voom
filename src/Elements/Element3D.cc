// -*- C++ -*-
//----------------------------------------------------------------------
//
//                           Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------


#include "Element3D.h"

namespace voom
{
   //! Constructor
  Element3D::Element3D( const NodeContainer & Nodes,
		        Material * Mat,
			Quadrature<3> & Quad,  
			Shape<3> & Sh, 
			double k): _nodes(Nodes), _k(k), _quad(&Quad), _sh(&Sh)
  { 
    for(ConstNodeIterator n = Nodes.begin(); n != Nodes.end(); n++) 
      _baseNodes.push_back(*n);
    
    //! initialize quad points
    _quadPoints.clear();
    
    for(Quadrature<3>::ConstPointIterator p = Quad.begin(); p!=Quad.end(); p++) {
      // compute shape functions at the Gauss point
      Sh.compute( p->coords );
      
      // create the QuadPointStruct
      _quadPoints.push_back( QuadPointStruct( p->weight, Mat, Sh, Nodes) ); 

      if( Sh.functions().size() != Nodes.size() ) {
	std::cout << "Number of nodes: " << Nodes.size() << std::endl
		  << "Number of functions: " << Sh.functions().size()
		  << std::endl
		  << "These should be equal." << std::endl;
	exit(0);
      }
    }

    _volume = 0.0;
    _strainEnergy = 0.0;
    compute(true, true, false);

    // compute initial element volume, which is given as sum of J*weight
    // however, quad point weight was already set as weight*J
    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++){
      _volume += p->weight;
    }      
    
  }


                            
  Element3D::QuadPointStruct::QuadPointStruct(double w,
					      Material * Mat,
				              Shape<3> & Sh,
					      const NodeContainer & Nodes): weight(w), material(Mat)
  {
    unsigned int a = 0, alpha = 0, i = 0;
    shapeFunctions = Sh.functions();

    const Shape<3>::DerivativeContainer & dnds = Sh.derivatives();
    
    // matrix jacobian dxds;
    Tensor3D dxds(0.0);
    for(i = 0; i < 3; i++) {
      for(alpha = 0; alpha < 3; alpha++) {	   
	for(a = 0; a < Nodes.size(); a++) {
	  dxds(i,alpha) +=  dnds[a](alpha)*( Nodes[a]->getPosition(i) );
	}
      }	      
    }
  
    // compute scalar jacobian and scale the quadrature weight with it
    // by calculating the determinant of dxds
    double jac = fabs(determinant(dxds));
   
    weight *= fabs(jac);
    // cout << weight << " jac  " << jac << endl;
    
    // invert matrix jacobian
    Tensor3D invJac(0.0);
    invert(dxds, invJac);
    
    // spatial derivatives dndx
    shapeDerivatives.resize(dnds.size());
    for(a = 0; a < Nodes.size(); a++) {	  
      for(i = 0; i < 3; i++) {
	shapeDerivatives[a](i) = 0.0;
	for(alpha = 0; alpha < 3; alpha++) {
	  shapeDerivatives[a](i) += dnds[a](alpha)*invJac(alpha,i);
	}
      }
    }
   
   
    return;
  }


  void Element3D::reset()
  {
    _volume = 0.0;
    
    QuadPointIterator p = _quadPoints.begin();
    Quadrature<3>::ConstPointIterator q = _quad->begin();

    for( ; p != _quadPoints.end(); p++, q++)
    {
       unsigned int a = 0, alpha = 0, i = 0;
       _sh->compute( q->coords );

       // reset shape functions at the quadrature point
       p->shapeFunctions = _sh->functions();

       const Shape<3>::DerivativeContainer & dnds = _sh->derivatives();
    
       // matrix jacobian dxds;
       Tensor3D dxds(0.0);
       for(i = 0; i < 3; i++) {
	 for(alpha = 0; alpha < 3; alpha++) {	   
	   for(a = 0; a < _nodes.size(); a++) {
	     dxds(i,alpha) +=  dnds[a](alpha)*( _nodes[a]->getPosition(i) );
	   }
	 }	      
       }
  
       // compute scalar jacobian and scale the quadrature weight with it
       // by calculating the determinant of dxds
       double jac = fabs(determinant(dxds));
   
       p->weight = q->weight*fabs(jac);
       _volume += p->weight;
       // cout << p->weight << " jac  " << jac << endl;
    
       // invert matrix jacobian
       Tensor3D invJac(0.0);
       invert(dxds, invJac);
    
       // spatial derivatives dndx
       for(a = 0; a < _nodes.size(); a++) {	  
	 for(i = 0; i < 3; i++) {
	   p->shapeDerivatives[a](i) = 0.0;
	   for(alpha = 0; alpha < 3; alpha++) {
	     p->shapeDerivatives[a](i) += dnds[a](alpha)*invJac(alpha,i);
	   }
	 }
       }

       
  
    }

    return;
  }


void Element3D::compute(bool f0, bool f1, bool f2)
  {
    // Initialize
    if( f0 ) {
      _energy = 0.0;
      _strainEnergy = 0.0;
    }

    blitz::Array< Vector3D, 1> _internalForce;
    if( f1 ) {
      _internalForce.resize( _nodes.size() );
      _internalForce = Vector3D(0.0);
    }

    unsigned int a = 0, i = 0, J = 0;
	
    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
    {  
      // compute deformation gradient
      Tensor3D F(0.0);
      
      const Shape<3>::DerivativeContainer &  DN = p->shapeDerivatives;

      for(a = 0; a < _nodes.size(); a++)
      {
	const Vector3D & xa = _nodes[a]->point();
	for(i = 0; i < 3 ; i++) {
	  for(J = 0; J < 3; J++) {
	    F(i,J) += xa(i)*DN[a](J);
	  } 
	}
      }

      if(determinant(F) > 0.0)
      {
	// send updated deformation gradient to material
	p->material->setDeformationGradient(F);
       
	// compute strain energy and/or 1st PK stress
	p->material->updateState(f0, f1, f2); 
      
	double weight = p->weight;

	// compute energy
	if ( f0 ){
	  // compute strain energy 
	  _strainEnergy += p->material->energyDensity()*weight;
	  if(_k > 0.0)
	    {
	      DeformationNode<3>::PositionVector X; 
	      DeformationNode<3>::Point x;
	      for(a = 0; a < _nodes.size(); a++)
		{
		  X = _nodes[a]->position();
		  x = _nodes[a]->point();
		  _strainEnergy += 0.5*_k*sqr(norm2(X-x))*weight;
		}
	    }
	}

	// compute forces
	if ( f1 ) 
	{
	    const Tensor3D & P = p->material->piolaStress();
	    // cout << P << endl;
	    // compute internal forces
	    // loop for all nodes to compute forces 
	    DeformationNode<3>::PositionVector X; 
	    DeformationNode<3>::Point x;
	    for (a = 0; a < _nodes.size(); a++) { 
	      X = _nodes[a]->position();
	      x = _nodes[a]->point();
	      for(i = 0; i < 3; i++) {
		for(J = 0; J < 3; J++) {
		  _internalForce(a)(i) += P(i,J)*DN[a](J);
		}  
		if (_k > 0.0)
		{  
		  _internalForce(a)(i) += _k*(x(i)-X(i)); 
		}
		_internalForce(a)(i) *= weight;
	      }
	    } // end nodes loop
	
	} // end force calcs
      
	// compute stiffness matrix
	if( f2 ) {
	  std::cerr << "No stiffness matrix " << std::endl;
	}

      } // end of if detF > 0 loop 
    
    }  // end quadrature loop 
    
    if(f0) {
      _energy = _strainEnergy;
    }

    // assemble element forces to nodes
    NodeIterator na;
    double f_ia = 0.0;
    if(f1) {
      a = 0;
      for(na = _nodes.begin();  na != _nodes.end(); na++, a++) {	
	for(i = 0; i < 3; i++) {
	  f_ia = _internalForce(a)(i);
	  (*na)->addForce( i, f_ia );
	}
      }
    }

    
    
  }



  vector<pair<Vector3D, vector<double > > >  Element3D::invariants(int &f)
  {
    vector<pair<Vector3D, vector<double > > > Invariants;
    unsigned int a = 0, i = 0, j =0, J = 0, k = 0, q = 0;

    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
    {  
      // compute deformation gradient
      Tensor3D F(0.0), C(0.0), Csquare(0.0);
      
      const Shape<3>::DerivativeContainer &  DN = p->shapeDerivatives;
      const Shape<3>::FunctionContainer &  N = p->shapeFunctions;
      Vector3D Location(0.0);

      for(a = 0; a < _nodes.size(); a++)
      {
	const Vector3D & xa = _nodes[a]->point();
	Location += N[a]*xa;

	for(i = 0; i < 3 ; i++) {
	  for(J = 0; J < 3; J++) {
	    F(i,J) += xa(i)*DN[a](J);
	  } 
	}
      }

      // Compute invariants
      double detF = determinant(F); 
      double detF_TwoThird = pow(detF, -2.0/3.0);

      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  for (k = 0; k < 3; k++) {
	    C(i,j) += F(k,i)*F(k,j);
	  }
	}
      }

      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  for (k = 0; k < 3; k++) {
	    Csquare(i,j) += C(i,k)*C(k,j);
	  }
	}
      }

      double trC = C(0,0) + C(1,1) + C(2,2);

      vector<double > InvLoc(3,0.0);
      InvLoc[0] = trC*detF_TwoThird;
      InvLoc[1] = 0.5*(pow(trC,2.0) - Csquare(0,0) - Csquare(1,1) - Csquare(2,2))*pow(detF_TwoThird, 2.0);
      InvLoc[2] = detF;//*detF; 
      q++;
      if(detF < 0.0)
      { 
	// std::cout<<"The jacobian is negative at ." << endl;
	f++;
	InvLoc[0] = -1.0;
	InvLoc[1] = -1.0;
	InvLoc[2] = -1.0;
      }

      
      
      Invariants.push_back(make_pair(Location, InvLoc));
      
    }

    return Invariants;
  }

} // namespace voom
