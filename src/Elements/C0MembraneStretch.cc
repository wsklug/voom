// -*- C++ -*-
//----------------------------------------------------------------------
//
//                      William S. Klug, Feng Feng
//                      Luigi Perotti
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

#include "C0MembraneStretch.h"
#include "ShapeTri3.h"

//#define _DEBUG_

namespace voom
{
  C0MembraneStretch::C0MembraneStretch(NodeContainer & Nodes,
				       ScalarFieldNode<3> * StretchNode,
				       ScalarFieldNode<3> * DirectionNode,
				       ScalarFieldNode<3> * phiNode,
				       const double &AngleOffset,
				       const Quadrature<2> & Quad,
				       EvansElastic_Stretch *Mat,
				       Shape<2> * shape,
				       bool InsertStretch,
				       bool InsertDirection,
				       bool InsertPhi)
    {
      _area = 0.0;
      _volume = 0.0;
      _energy = 0.0;
      _confEnergy = 0.0;
      //! initialize nodes
      _baseNodes.insert(_baseNodes.begin(), Nodes.begin(), Nodes.end());
      if(InsertStretch) {_baseNodes.push_back(StretchNode);};
      if(InsertDirection) {_baseNodes.push_back(DirectionNode);};
      if(InsertPhi) {_baseNodes.push_back(phiNode);};
      _nodes = Nodes;
      
      _stretchNode = StretchNode;
      _directionNode = DirectionNode;
      _phiNode = phiNode;
      _angleOffset = AngleOffset;
	  
      //! initialize materials and shape functions
      for( Quadrature<2>::ConstPointIterator p = Quad.begin();  p != Quad.end(); p++)
      {
	shape->compute(p->coords);
 	const Shape<2>::FunctionContainer   & fun = shape->functions();
        const Shape<2>::DerivativeContainer & der = shape->derivatives();
	_quadPoints.push_back( QuadPointStruct(p->weight, Mat, fun, der) ); 
       
	assert( fun.size() == Nodes.size() );
      }

      // Update reference configuration and initiliaze values
      updateRefConfiguration(); 

      // cout <<  _baseNodes.size() << endl;

    }
   


  void C0MembraneStretch::compute(bool f0, bool f1, bool f2)
  {
    
    int b = 0, alpha = 0, i = 0, ib = 0;
    Vector3D zero(0);
    //! forces on the element - temp variable
    blitz::Array< Vector3D, 1> _internalForce;
    _internalForce.resize(_nodes.size() );

    if( f0 ) {
      _energy = 0.0;
      _confEnergy = 0.0;
    }

    if( f1 ) {
      _internalForce = zero;
    }
    
    // loop for every quadrature point
    for(QuadPointIterator p = _quadPoints.begin(); p!=_quadPoints.end(); p++)
    {
      // compute Shell geometry
      		
      // compute position, basis vector, and derivatives of basis vector
      Vector3D x(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      a = zero, zero;
      aPartials = zero, zero, zero, zero;

      const Shape<2>::FunctionContainer & N = p->shapeFunctions;
      const Shape<2>::DerivativeContainer & DN = p->shapeDerivatives;

      for (b = 0; b < _nodes.size(); b++)
      {
	const DeformationNode<3>::Point & xb = _nodes[b]->point();
	x 	       +=   N[b]    * xb;
	a(0)           +=  DN[b](0) * xb;
	a(1)           +=  DN[b](1) * xb;
      }	
      
      // compute shell geometry
      ShellGeometry geometry(a , aPartials);
      const Vector3D& d = geometry.d();
      const tvmet::Vector< Vector3D, 2 >& aDual = geometry.aDual();
      		
      // material part
      EvansElastic_Stretch *material = p->material;
      // store the deformed geometry in shell geometry class
      material->setGeometry(geometry);
      material->setEta(_stretchNode->point());
      material->setTheta(_angleOffset + _directionNode->point() );
      if (_phiNode != NULL) {
	material->setPhi(_phiNode->point() );
      }

      // compute strain energy, stress and moment resultants
      material->updateState(f0, f1, f2); 
      		
      // 0.5 is not needed
      const double metric = geometry.metric();
      const double refMetric = ( material->refShellGeometry()).metric();
      const double jacobian = metric/refMetric;
      const double weight =  metric * p->weight;

      // compute area for the area constaint energy
      _area += weight;

      // compute volume for the volume constraint energy
      _volume +=  dot(d,x) * weight / 3.0;
      
      // compute energy
      if ( f0 ){
	// compute strain energy 
	_energy += material->energyDensity() * weight;
	_confEnergy += material->conformationalEnergy() * weight;
      }
      
      // compute forces
      if ( f1 ) {
	// get stress and moment resultants
	const tvmet::Vector< Vector3D, 3 >& sr = material->stressResultants();
	const double stretchForce = material->stretchForce();
	const double directionForce = material->directionForce();
	const double PhiForce = material->phiForce();
	
	// loop for all nodes to compute forces 
	for (b = 0; b < _nodes.size(); b++) 
	{
	  // compute internal forces
	  // calculate the gradient of the derivatives of the director w.r.t curvilinear coords
	  for (alpha = 0; alpha < 2; alpha++)
	  {
	    // Stress Resultant part
	    Vector3D ftmp;
	    ftmp = sr(alpha) *  DN[b](alpha) * weight;
 	    _internalForce(b)  += ftmp;
	  }
	} // end nodes loop
	
	_stretchNode->addForce(stretchForce*weight);
	_directionNode->addForce(directionForce*weight);
	if (_phiNode != NULL) {
	  _phiNode->addForce(PhiForce*weight);
	}
        // cout << "sheraNode->force() = " << _stretchNode->force() << endl;
	
      } // end force calcs

      // compute stiffness matrix
      if( f2 ) {//assert(0)
	       } // Stiffness is not available }
      
    } // end quadrature loop

    if( f1 ) {
      b = 0;
      for(NodeIterator nb = _nodes.begin();  nb != _nodes.end(); nb++, b++)
	for(i = 0; i < 3; i++)
	{
	    (*nb)->addForce( i, _internalForce(b)(i) );
	}
    }
   

  }



  void C0MembraneStretch::invariants(double& I1, double& J)
  {
    int b = 0, i = 0, j = 0, k = 0, alpha = 0;
    Vector3D zero(0);
    
    // Loop for every quadrature point
    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
    {
      // compute position, basis vector, and derivatives of basis vector
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      a = zero, zero;
      aPartials = zero, zero, zero, zero;

      const Shape<2>::FunctionContainer   &  N = p->shapeFunctions;
      const Shape<2>::DerivativeContainer & DN = p->shapeDerivatives;

      for (b = 0; b < _nodes.size(); b++)
      {
	const DeformationNode<3>::Point & xb = _nodes[b]->point();
	a(0) += DN[b](0) * xb;
	a(1) += DN[b](1) * xb;
      }

      // compute shell geometry
      ShellGeometry defgeometry(a, aPartials);
      a = zero, zero;
      aPartials = zero, zero, zero, zero;
      for (b = 0; b < _nodes.size(); b++)
      {
	const DeformationNode<3>::Point & Xb = _nodes[b]->position();
	a(0) +=  DN[b](0) * Xb;
	a(1) +=  DN[b](1) * Xb;
      }
      ShellGeometry refgeometry( a, aPartials );
      
      const tvmet::Vector<Vector3D, 2> & basis   = defgeometry.a();
      const tvmet::Vector<Vector3D, 2> & refDual = refgeometry.aDual();

      Tensor3D F(0.0);
      for(alpha = 0; alpha<2; alpha++){
        for(i=0; i<3; i++) {
  	  for(j=0; j<3; j++) {
	    F(i,j) += basis(alpha)(i)*refDual(alpha)(j);	  
	  }
        }
      }

      Tensor3D C(0.0);
      C = trans(F)*F;
      I1 = C(0,0)+C(1,1)+C(2,2);
      double trCSquare = 0.0;
      for(i = 0; i < 3; i++) {
        for(k = 0; k < 3; k++) {
	  trCSquare += C(i,k)*C(k,i);
        }
      }
      J = sqrt( 0.5*(I1*I1 - trCSquare) );
    }
   

  }



  Vector3D C0MembraneStretch::PushForwardOperator(Vector3D & Nbar)
  {  
    int b = 0, i = 0, j = 0, alpha = 0;
    Vector3D zero(0), nbar(0);
    
    // Loop over quadrature point
    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
    {
      // compute position, basis vector, and derivatives of basis vector
      tvmet::Vector< Vector3D, 2 > a;
      a = zero, zero;

      const Shape<2>::DerivativeContainer & DN = p->shapeDerivatives;

      for (b = 0; b < _nodes.size(); b++)
      {
	const DeformationNode<3>::Point & xb = _nodes[b]->point();
	a(0) += DN[b](0) * xb;
	a(1) += DN[b](1) * xb;
      }

      const tvmet::Vector<Vector3D, 2> & refDual = p->material->refShellGeometry().aDual();

      Tensor3D F(0.0);
      for(alpha = 0; alpha<2; alpha++){
        for(i=0; i<3; i++) {
  	  for(j=0; j<3; j++) {
	    F(i,j) += a(alpha)(i)*refDual(alpha)(j);	  
	  }
        }
      }

      // nbar += F*Nbar*(p->weight); 
      nbar = F*Nbar;
    }
    
    nbar /= norm2(nbar);
   
    return nbar;
  }



  void C0MembraneStretch::updateRefConfiguration() 
  {
    int b = 0;
    Vector3D zero(0);

    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++)
    {
      // compute Shell geometry			
      // compute position, basis vector, and derivatives of basis vector
      Vector3D X(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      a = zero, zero;
      aPartials = zero, zero, zero, zero;
			
      const Shape<2>::FunctionContainer   &  N = p->shapeFunctions;
      const Shape<2>::DerivativeContainer & DN = p->shapeDerivatives;

      for (b = 0; b < _nodes.size(); b++)
      {
	const DeformationNode<3>::Point & Xb = _nodes[b]->position();
	X    +=  N[b]    * Xb;
	a(0) += DN[b](0) * Xb;
	a(1) += DN[b](1) * Xb;
      }
			
      // compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
			
      // store the reference geometry in shell geometry class
      p->material->setRefGeometry(refgeometry);
    }    
    compute(true,true,false);
   
  }
   


  // added to set the reference configuration explicitly
  void C0MembraneStretch::SetRefConfiguration(double edgelen) 
  {
    for(QuadPointIterator p = _quadPoints.begin(); p!=_quadPoints.end(); p++)
    {
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      a(0) = edgelen, 0.0, 0.0;
      a(1) = edgelen*cos(M_PI/3.0), edgelen*sin(M_PI/3.0), 0.0;

      // compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
			
      // store the reference geometry in shell geometry class
      p->material->setRefGeometry(refgeometry);

    }    
    compute(true,true,false);
    
   
  }



  Vector3D C0MembraneStretch::computePosition(const double s1, const double s2)
  {
    // CornerValences v(6,6,6);
    tvmet::Vector<unsigned int, 3> v(6,6,6);
    // Array1D paraCoords(2);
    ShapeTri3::CoordinateArray parCoord;
    parCoord = s1, s2;
    const int nodes = v(0) + v(1) + v(2) - 6;
    //
    // create a loop shell shape function object
    ShapeTri3 shp( parCoord );
    //
    // new position
    Vector3D pos(0.0);
    NodeIterator p = _nodes.begin();
    for(; p != _nodes.end(); p ++){
      const int i = std::distance(_nodes.begin(), p);
      pos += shp.functions()[i] * (*p)->point();
    }
		
    return pos;
  }


  const Tensor3D C0MembraneStretch::cauchyStress()
  {
    Tensor3D sigma(0.0);
    QuadPointIterator QPit = _quadPoints.begin();
    double Atot = 0.0;
    for( ; QPit != _quadPoints.end(); QPit++)
    {
      sigma += QPit->material->cauchyStress()*QPit->weight;
      Atot += QPit->weight;
    }
    sigma /= Atot;
    // cout << sigma << " " << Atot << endl;
    
    return sigma;
  };
  
  const vector<double > C0MembraneStretch::matInvariants()
  {
    vector<double > invariants(2, 0.0), tempInv(2, 0.0);
    QuadPointIterator QPit = _quadPoints.begin();
    double Atot = 0.0;
    for( ; QPit != _quadPoints.end(); QPit++)
    {
      tempInv = QPit->material->invariants();
      invariants[0] += tempInv[0]*QPit->weight;
      invariants[1] += tempInv[1]*QPit->weight;
      Atot += QPit->weight;
    }
    invariants[0] /= Atot;
    invariants[1] /= Atot;

    return invariants;
  };
  
} // namespace voom
