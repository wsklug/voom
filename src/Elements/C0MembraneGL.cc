// -*- C++ -*-
//----------------------------------------------------------------------
//
//                      William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

#include "C0MembraneGL.h"

////////////////////////////////////////////////////////////////////
// LAPACK subroutine for computing eigenvalues and eigenvectors
//
// define prototype of LAPACK functions
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		       double *w, double *work, int *lwork, int *info);

//#define _DEBUG_

namespace voom
{

  C0MembraneGL::QuadPointStruct::QuadPointStruct
  (double w, const MaterialType * m, Shape<2> * s) 
    : weight(w), 
      material(*m), 
      shapeFunctions( s->functions() ), 
      shapeDerivatives( s->derivatives() ) 
  {}

  C0MembraneGL::C0MembraneGL
  ( const DefNodeContainer  & defNodes,
    const GLNodeContainer   & glNodes ,
    const MaterialType * material,
    Quadrature<2> * quad,
    Shape<2> * shape )
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
    _volume = 0.0;
    _area = 0.0;

    //! initialize materials and shape functions
    _quadPoints.clear();
    for(Quadrature<2>::ConstPointIterator p=quad->begin(); 
	p!=quad->end(); p++) {
      // compute shape functions at the Gauss point
      shape->compute( p->coords );
      
      // create the QuadPointStruct
      _quadPoints.push_back( QuadPointStruct( p->weight, material, shape ) ); 
    }

    updateRefConfiguration(); 

    //! initialize values = 0 by default
    //compute(false, false, false); 
    compute(true, true, false); 
    
    return;
    
  }

  void C0MembraneGL::updateRefConfiguration() {
    for(QuadPointIterator p=_quadPoints.begin(); 
	p!=_quadPoints.end(); p++){
      //
      // compute Shell geometry
      //
			
      // compute position, basis vector, and derivatives of basis vector
      tvmet::Vector< Vector3D, 2 > a;
      tvmet::Matrix< Vector3D, 2, 2 > aPartials;
      Vector3D zero(0.0);
      a = zero, zero;
      aPartials = 
	zero, zero, 
	zero, zero;
			
      const Shape<2>::FunctionContainer & N  = p->shapeFunctions;
      const Shape<2>::DerivativeContainer & DN = p->shapeDerivatives;

      for (int b = 0; b < _defNodes.size(); b++){

	const Vector3D & Xb = _defNodes[b]->position();
	a(0)           +=  DN[b](0)   * Xb;
	a(1)           +=  DN[b](1)   * Xb;
      }

      // compute shell geometry
      ShellGeometry refgeometry( a, aPartials );
			
      // store the reference geometry in shell geometry class
      p->material.setRefGeometry(refgeometry);

      //material.updateState(true, true, true);

    }    
    compute(true,true,true);
    
    return;
  }
   


  void C0MembraneGL::compute(bool f0, bool f1, bool f2)
  {
    //
    // initialize things
    //
    // if applied volume constraint, _volume is needed for computing
    // if applied area constraint, _area is needed for computing
    // both energy and forces
    _volume = 0.0;
    _area = 0.0;
    // _trC = 0.0;

   
		
    if( f0 ) {
      _energy = 0.0;
      _strainEnergy = 0.0;
      _work = 0.0;
    }

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
	const DeformationNode<3>::Point & xb = _defNodes[b]->point();
	x 	       +=   N[b]     * xb;
	a(0)           +=  DN[b](0)   * xb;
	a(1)           +=  DN[b](1)   * xb;
      }
	
      // compute GL field and gradient

      double c=0.0;
      Vector2D dc( 0.0 );

      for (int b = 0; b < _glNodes.size(); b++){
	double cb = _glNodes[b]->point();
	c 	       +=   N[b]      * cb;
	dc(0)           +=  DN[b](0)   * cb;
	dc(1)           +=  DN[b](1)   * cb;
      }

      if(norm2(a(0)) > 1.0e3) {
	std::cout << "x = " << x << std::endl
	     << "a(0) = " << a(0) << std::endl
	     << "a(1) = " << a(1) << std::endl;
	
	for (int b = 0; b < _defNodes.size(); b++){
	  const DeformationNode<3>::Point & xb = _defNodes[b]->point();
	  std::cout << "x_" << b << " = " << xb << std::endl;
	}
	exit(0);
      }
	
      // compute shell geometry
      ShellGeometry geometry( a, aPartials );
      const Vector3D& d = geometry.d();
      const tvmet::Vector< Vector3D, 2 >& aDual = geometry.aDual();
			
      // store the deformed geometry in shell geometry class
      p->material.setGeometry(geometry);

      p->material.setField( c, dc );

      // compute strain energy, stress and moment resultants
      p->material.updateState(f0, f1, f2); 
			
      const double metric =  geometry.metric();
      const double refMetric = ( p->material.refShellGeometry()).metric();
      const double jacobian = metric/refMetric;
      const double weight =  metric * p->weight;

      // compute area for the area constaint energy
      _area += weight;

      // compute volume for the volume constraint energy
      _volume +=  dot(d,x) * weight / 3.0;
      
      // compute energy
      if ( f0 ){

	// compute strain energy 
	_strainEnergy += p->material.energyDensity() * weight;
              
      }

      // compute forces
      if ( f1 ) {
	//
	// get stress and moment resultants
	const tvmet::Vector< Vector3D, 3 >& n = p->material.stressResultants();

	const double mu = p->material.GLResultants();
	const tvmet::Vector< double, 2 >& lambda = p->material.GLGradientResultants();

	// loop for all nodes to compute forces 
	for (int a=0; a<_defNodes.size(); a++) {

	  // compute internal forces

	  // calculate the gradient of the derivatives of the director
	  // w.r.t curvilinear coords
	  for ( int alpha = 0; alpha < 2; alpha++){
	    // Stress Resultant part
	    Vector3D f;
	    f = n(alpha) *  DN[a](alpha) * weight;
	    for(int i=0; i<3; i++) _defNodes[a]->addForce( i, f(i) );
	  } 

	  // nodal chemical potential
	  _glNodes[a]->addForce( mu*N[a]*weight );
	  for(int beta=0; beta<2; beta++){
	    _glNodes[a]->addForce( lambda(beta)*DN[a](beta)*weight );
	  } 

	} // end nodes loop

	
      } // end force calcs

    } // end quadrature loop


    if(f0) {
      _energy = _strainEnergy; //- _work + tension * _area;
    }

  }

  
} // namespace voom
