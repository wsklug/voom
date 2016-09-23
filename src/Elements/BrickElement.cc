// -*- C++ -*-
//----------------------------------------------------------------------
//
//                           Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#include "./BrickElement.h"

////////////////////////////////////////////////////////////////////
// LAPACK subroutine for computing eigenvalues and eigenvectors
//
// define prototype of LAPACK functions
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
		       double *w, double *work, int *lwork, int *info);

//#define _DEBUG_

namespace voom
{
                                
  BrickElement::QuadPointStruct::QuadPointStruct(double w,
						 MaterialType * m,
						 Shape<3> * s,
						 const NodeContainer & nds): weight(w), material(m)
  {
    unsigned int a = 0, alpha = 0, i = 0;
    shapeFunctions = s->functions();

    const Shape<3>::DerivativeContainer & dnds = s->derivatives();
 
    // matrix jacobian dxds;
    Tensor3D dxds(0.0);
     for(i = 0; i < 3; i++) {
      for(alpha = 0; alpha < 3; alpha++) {	   
	for(a = 0; a < nds.size(); a++) {
	  dxds(i,alpha) +=  dnds[a](alpha)*( nds[a]->getPosition(i) );
	}
      }	      
    }

    // compute scalar jacobian and scale the quadrature weight with it
    // by calculating the determinant of dxds
    double jac = determinant(dxds);
   
    weight *= jac;
    
    // invert matrix jacobian
    Tensor3D invJac(0.0);
    invert(dxds, invJac);
    
    // spatial derivatives dndx
    shapeDerivatives.resize(dnds.size());
    for(a = 0; a < nds.size(); a++) {	  
      for(i = 0; i < 3; i++) {
	shapeDerivatives[a](i) = 0.0;
	for(alpha = 0; alpha < 3; alpha++) {
	  shapeDerivatives[a](i) += dnds[a](alpha)*invJac(alpha,i);
	}
      }
    }
    
    return;
  }


  //! Constructor
  BrickElement::BrickElement( const NodeContainer & nodes,
			      MaterialType * mat,
			      Quadrature<3> * quad,  
			      Shape<3> * shape )
  {
    //! initialize NodeContainer
    unsigned nNodes = nodes.size();
      
    _nodes = nodes;
    
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
      _baseNodes.push_back(*n);
    
    //! initialize quad points
    _quadPoints.clear();

    for(Quadrature<3>::ConstPointIterator p=quad->begin(); p!=quad->end(); p++) {
      // compute shape functions at the Gauss point
      shape->compute( p->coords );
      
      // create the QuadPointStruct
      _quadPoints.push_back( QuadPointStruct( p->weight, mat, shape, _nodes) ); 
      if( shape->functions().size() != nNodes ) {
	std::cout << "Number of nodes: " << nNodes << std::endl
		  << "Number of functions: " << shape->functions().size()
		  << std::endl
		  << "These should be equal." << std::endl;
	exit(0);
      }
    }

    //! allocate memory for mechanics variables
    _internalForce.resize( nNodes );
    
    //! initialize values = 0 by default
    compute(true, true, false); 

    _volume = 0.0;

    // compute initial element volume, which is given as sum of J*weight
    // however, quad point weight was already set as weight*J
    for(QuadPointIterator p = _quadPoints.begin(); p != _quadPoints.end(); p++){
      _volume += p->weight;
    }      
    
  }
 
void BrickElement::compute(bool f0, bool f1, bool f2)
  {
    // Initialize
    if( f0 ) {
      _energy = 0.0;
      _strainEnergy = 0.0;
    }

    if( f1 ) {
      _internalForce = Vector3D(0.0);
    }

    unsigned int a = 0, i = 0, J = 0;

    for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++)
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
      
      MaterialType * material = p->material;
      
      // send updated deformation gradient to material
      material->setDeformationGradient(F);
      
      // compute strain energy and/or 1st PK stress
      material->updateState(f0, f1, f2); 
      
      double weight = p->weight;

      // compute energy
      if ( f0 ){
	// compute strain energy 
	_strainEnergy += material->energyDensity()*weight;
      }
      
      // compute forces
      if ( f1 ) 
      {
	const Tensor3D & P = material->piolaStress();
	// compute internal forces
	// loop for all nodes to compute forces 
	for (a = 0; a < _nodes.size(); a++) {
	  for(i = 0; i < 3; i++) {
	    for(J = 0; J < 3; J++) {
	      _internalForce(a)(i) += P(i,J)*DN[a](J);
	    }  
	    _internalForce(a)(i) *= weight;
	  }
	} // end nodes loop
	
      } // end force calcs
      
	// compute stiffness matrix
      if( f2 ) {
        std::cerr << "No stiffness matrix " << std::endl;
      }

    } // end quadrature loop
    
    
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



} // namespace voom
