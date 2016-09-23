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

namespace voom
{

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
      
      Material & material = p->material;
      
      // send updated deformation gradient to material
      material.setDeformationGradient(F);
      
      // compute strain energy and/or 1st PK stress
      material.updateState(f0, f1, f2); 
      
      double weight = p->weight;
           
      // compute energy
      if ( f0 ){
	// compute strain energy 
	_strainEnergy += material.energyDensity()*weight;
      }
      
      // compute forces
      if ( f1 ) 
      {
	const Tensor3D & P = material.piolaStress();
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





  double BrickElement::CalcMisesStress()
  {	
    double vMStress = 0.0;
    unsigned int numQP = 0, a = 0, i = 0, j = 0, J = 0, K = 0;

    for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++, numQP++)
    {
      // compute deformation gradient
      Tensor3D F(0.0);
      
      const Shape<3>::DerivativeContainer &  DN = p->shapeDerivatives;
      
      for(a = 0; a < _nodes.size(); a++)
      {
	const Vector3D & xa = _nodes[a]->point();
	for(i = 0; i < 3; i++) {
	  for(J = 0; J < 3; J++) {
	    F(i,J) += xa(i)*DN[a](J);
	  } 
	}
      }
      
      Material & material = p->material;
      
      // send updated deformation gradient to material
      material.setDeformationGradient(F);
      
      // compute 1st PK stress
      material.updateState(false, true, false); 
      
      const Tensor3D & P = material.piolaStress();

      double jac = determinant(F);

      Tensor3D Cauchy(0.0);
      // = 1.0/jac*P*transF;
      // Tensor3D  Cauchy = P*(tvmet::trans(F)); // (P*(tvmet::trans(F)))/jac;

      /*
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  for(K = 0; K < 3; K++) {
	    Cauchy(i,j) += 1.0/jac*P(i,K)*F(j,K);
	  }
	}
      }
      */
      vMStress += 0.707106781186548*sqrt((Cauchy(0,0)-Cauchy(1,1))*(Cauchy(0,0)-Cauchy(1,1)) + 
					 (Cauchy(1,1)-Cauchy(2,2))*(Cauchy(1,1)-Cauchy(2,2)) + 
					 (Cauchy(2,2)-Cauchy(0,0))*(Cauchy(2,2)-Cauchy(0,0)) + 
					 6.0*(Cauchy(0,1)*Cauchy(0,1)+Cauchy(0,2)*Cauchy(0,2)+Cauchy(1,2)*Cauchy(1,2)));
      
    } // end quadrature loop
    
    vMStress /= double(numQP);
    
    return vMStress;    
  }




  

  Tensor3D BrickElement::CalcCauchyStress()
  {	
    unsigned int numQP = 0, a = 0, i = 0, j = 0, J = 0, K = 0;

    Tensor3D Cauchy(0.0);
    
    for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++, numQP++)
    {
      // compute deformation gradient
      Tensor3D F(0.0);
      
      const Shape<3>::DerivativeContainer &  DN = p->shapeDerivatives;
      
      for(a = 0; a < _nodes.size(); a++)
      {
	const Vector3D & xa = _nodes[a]->point();
	for(i = 0; i < 3; i++) {
	  for(J = 0; J < 3; J++){
	    F(i,J) += xa(i)*DN[a](J);
	  } 
	}
      }
      
      Material & material = p->material;
      
      // send updated deformation gradient to material
      material.setDeformationGradient(F);
      
      // compute strain energy and/or 1st PK stress
      material.updateState(false, true, false); 
      
      const Tensor3D & P = material.piolaStress();

      double jac = determinant(F);

      // Tensor3D Cauchy(0.0);
      // = 1.0/jac*P*transF;
      Cauchy = (P*(tvmet::trans(F)))/jac;
      
      /*
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  for(K = 0; K < 3; K++) {
	    Cauchy(i,j) += 1.0/jac*P(i,K)*F(j,K);
	  }
	}
      }
      */
    } // end quadrature loop
    
    Cauchy /= double(numQP);

    return Cauchy;    
  }






  Vector3D BrickElement::CalcPrincipalStrains()
  {	

    Vector3D prinStrain;
    int numQP = 0, a = 0, i = 0, J = 0;
    
    for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++, numQP++)
    {
      // compute deformation gradient
      Tensor3D F(0.0);
      
      const Shape<3>::DerivativeContainer &  DN = p->shapeDerivatives;
      
      for(a = 0; a < _nodes.size(); a++) {
	const Vector3D & xa = _nodes[a]->point();
	for(i = 0; i < 3; i++) {
	  for(J = 0; J < 3; J++) {
	    F(i,J) += xa(i)*DN[a](J);
	  } 
	}
      }
      
      Tensor3D C(3,3);
      // C:  right Cauchy-green strain tensor
      // Tensor3D C = (tvmet::trans(F))*F;

      // compute Eigenvalues and Eigenvectors by calling LAPACK library
      char jobz = 'N';
      char uplo = 'L';
      int  n    = 3;
      int  lda  = n;
      int  lwork = 3*n-1;
      int  info;
      double eigenvalues[3];
      double work[lwork];
      
      // calling lapack function here to compute
      // eigenvalues of C1
      dsyev_(&jobz, &uplo, &n, C.data(),
	     &lda, eigenvalues, work, &lwork, &info);

      for(int i=0; i<3; i++) {
	prinStrain(i) = eigenvalues[i];
      }

      if (info != 0) {
	std::cout << "Something is wrong in DSYEV_" << std::endl;
	exit(0);
      }
      
    } // end quadrature loop
    
    prinStrain /= double(numQP);

    return prinStrain;    
  }






  BrickElement::QuadPointStruct::
  QuadPointStruct(double w,
		  const MooneyRivlin & m,
		  const Shape<3> & s, 
		  const BrickElement::NodeContainer & nds): weight(w), material(m)
  {
    unsigned int a = 0, alpha = 0, i = 0;
    shapeFunctions = s.functions();
    
    // Compute spatial derivatives of shape functions from
    // parametric derivatives by transforming with matrix jacobian
    // of isoparametric mapping.
    
    // parametric derivatives from shape object
    const Shape<3>::DerivativeContainer & dnds = s.derivatives();
    
    // matrix jacobian dxds;
    Tensor3D dxds(0.0);
     for(i = 0; i < 3; i++) {
      for(alpha = 0; alpha < 3; alpha++) {	   
	for(a = 0; a < nds.size(); a++) {
	  dxds(i,alpha) += dnds[a](alpha)*( nds[a]->getPosition(i) );
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
  BrickElement::BrickElement( const Quadrature<3> & quad,
			      const Material & mat,
			      const NodeContainer & nodes,
			      const Shape<3> & shape
			      )
  {
    //! initialize NodeContainer
    unsigned nNodes = nodes.size();
      
    _nodes = nodes;
    
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
      _baseNodes.push_back(*n);
    
    //! initialize quad points
    _quadPoints.clear();

    std::vector<Point> Points = quad.points();

    for(unsigne int p=0; p<Coords.size(); p++) {
      shape->compute(Points[p].coords);
      _quadPoints.push_back( QuadPointStruct(Points[p].weight, mat, shpape, _nodes) ); 
      if( shp.functions().size() != nNodes ) {
	std::cout << "Number of nodes: " << nNodes << std::endl
		  << "Number of functions: " << shp.functions().size()
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

} // namespace voom
