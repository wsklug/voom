// -*- C++ -*-
//----------------------------------------------------------------------
//
//                      William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

#include "MembraneGLImplicitMass.h"

////////////////////////////////////////////////////////////////////
// LAPACK subroutine for computing eigenvalues and eigenvectors
//
// define prototype of LAPACK functions
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		       double *w, double *work, int *lwork, int *info);

//#define _DEBUG_

namespace voom
{

  MembraneGLImplicitMass::QuadPointStruct::QuadPointStruct
  (double w, Shape<2> * s) 
    : weight(w), 
      shapeFunctions( s->functions() ), 
      shapeDerivatives( s->derivatives() ) 
  {}

  MembraneGLImplicitMass::MembraneGLImplicitMass
  ( const DefNodeContainer  & defNodes,
    const GLNodeContainer   & glNodes ,
    Quadrature<2> * quad,
    Shape<2> * shape,
    double Mx0, double Mx1, double Meta0 )
  {

    
    if( shape->functions().size() != defNodes.size() ) {
      std::cout << "Number of deformation nodes: " << defNodes.size() 
		<< std::endl
		<< "Number of functions: " << shape->functions().size()
		<< std::endl
		<< "These should be equal." << std::endl;
      exit(0);
    }

    if( shape->functions().size() != glNodes.size() ) {
      std::cout << "Number of GL nodes: " << glNodes.size() 
		<< std::endl
		<< "Number of functions: " << shape->functions().size()
		<< std::endl
		<< "These should be equal." << std::endl;
      exit(0);
    }
    
    
    _Mx0 = Mx0;
    _Mx1 = Mx1;
    _Meta0 = Meta0;

    //! initialize NodeContainers
    _defNodes = defNodes;
    _glNodes = glNodes;
    
    for(ConstDefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
      assert( *n != 0 );
      _baseNodes.push_back(*n);
    }    
    for(ConstGLNodeIterator n=_glNodes.begin(); n!=_glNodes.end(); n++) {
      assert( *n != 0 );
      _baseNodes.push_back(*n);
    }    

    int nNodes = _defNodes.size();

    _Mx.resize( nNodes, nNodes);
    _Meta.resize( nNodes, nNodes);

    _Mx = 0.0;
    _Meta = 0.0;

    _x_prev.resize(nNodes);

    _eta_prev.resize(nNodes);
    //! initialize materials and shape functions
    _quadPoints.clear();
    for(Quadrature<2>::ConstPointIterator p=quad->begin(); 
	p!=quad->end(); p++) {
      // compute shape functions at the Gauss point
      shape->compute( p->coords );
      
      // create the QuadPointStruct
      _quadPoints.push_back( QuadPointStruct( p->weight, shape ) ); 
    }

    //! initialize values = 0 by default
    //compute(false, false, false); 

    return;
    
  }

   


  void MembraneGLImplicitMass::step(double dt)
  {
    // recompute masses

    _Mx = 0.0;
    _Meta = 0.0;

    // loop for every quadrature point
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
      //
      // compute Shell geometry
      //
			
      // compute position, basis vector, and derivatives of basis vector
      Vector3D x(0.0);
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      Vector3D zero(0);
      a = zero, zero;
      aPartials = zero, zero, zero, zero;

      const Shape<2>::FunctionContainer & N 
	= p->shapeFunctions;
      const Shape<2>::DerivativeContainer & DN 
	= p->shapeDerivatives;

      for (int b = 0; b < _defNodes.size(); b++){
	const DeformationNode<3>::Point & xb = _defNodes[b]->position();
	x 	       +=   N[b]     * xb;
	a(0)           +=  DN[b](0)   * xb;
	a(1)           +=  DN[b](1)   * xb;
      }
	
      // compute GL field and gradient

      double eta=0.0;

      for (int b = 0; b < _glNodes.size(); b++){
	double eta_b = _glNodes[b]->point();
	eta += N[b] * eta_b;
      }

      if(norm2(a(0)) > 1.0e3) {
	std::cout << "x = " << x << std::endl
	     << "a(0) = " << a(0) << std::endl
	     << "a(1) = " << a(1) << std::endl;
	
	for (int b = 0; b < _defNodes.size(); b++){
	  const DeformationNode<3>::Point & xb = _defNodes[b]->position();
	  std::cout << "x_" << b << " = " << xb << std::endl;
	}
	exit(0);
      }
	
      // compute shell geometry
      ShellGeometry geometry( a, aPartials );
      const Vector3D& d = geometry.d();
      const tvmet::Vector< Vector3D, 2 >& aDual = geometry.aDual();
			
      const double metric =  geometry.metric();
      const double weight =  metric * p->weight;
      
      double Mx = _Mx0 + _Mx1*eta;
      // loop for all nodes to compute forces 
      for (int a=0; a<_defNodes.size(); a++) {
	for (int b=a; b<_defNodes.size(); b++) {
	  // position mass
	  _Mx(a,b) += Mx*N[a]*N[b]*weight;
	  //symetry
	  _Mx(b,a) = _Mx(a,b);

	  // GL mass
	  _Meta(a,b) += _Meta0*N[a]*N[b]*weight;
	  //symetry
	  _Meta(b,a) = _Meta(a,b);
	}
      } // end nodes loop
      

    } // end quadrature loop

    _Mx /= dt;
    _Meta /= dt;

    // save positions and order parameters
    for(int a=0; a<_defNodes.size(); a++ ) {  
      _x_prev[a] = _defNodes[a]->point();
      _eta_prev[a] = _glNodes[a]->point();
    }

  }

  void MembraneGLImplicitMass::compute(bool f0, bool f1, bool f2)
  {
    if( f0 ) {
      _energy = 0.0;
      _strainEnergy = 0.0;
      
      for (int a=0; a<_defNodes.size(); a++) {

	Vector3D dx_a(0.0);
	dx_a = _defNodes[a]->point() - _x_prev[a];

	double deta_a = _glNodes[a]->point() - _eta_prev[a];

	for (int b=0; b<_defNodes.size(); b++) {

	  Vector3D dx_b(0.0);
	  dx_b = _defNodes[b]->point() - _x_prev[b];

	  double deta_b = _glNodes[b]->point() - _eta_prev[b];

	  _energy += 0.5*_Mx(a,b)*dot(dx_a,dx_b) + 0.5*_Meta(a,b)*deta_a*deta_b;
	}
      }
    }
    
    // compute forces
    if ( f1 ) {
      for (int a=0; a<_defNodes.size(); a++) {

	Vector3D f_a(0.0);
	double mu_a(0.0);
	for (int b=0; b<_defNodes.size(); b++) {
	  f_a += _Mx(a,b)*(_defNodes[b]->point() - _x_prev[b]);
	  mu_a += _Meta(a,b)*(_glNodes[b]->point() - _eta_prev[b]);
	}

	for(int i=0; i<3; i++) {
	  _defNodes[a]->addForce(i,f_a(i));
	}
	// nodal chemical potential
	_glNodes[a]->addForce( mu_a );
	
      } // end nodes loop
      
    } // end force calcs


  }

  
} // namespace voom
