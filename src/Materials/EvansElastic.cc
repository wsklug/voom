// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/09/28 03:50:48  ma
// Adding cytoskeleton part to SCElasttic, using Evans model
// to claculate cytoskeleton energy
//
// Revision 1.9  2005/06/18 17:05:02  klug
// Added Gaussian curvature term to energy.
//
// Revision 1.8  2005/05/23 17:38:39  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file EvansElastic.cc

  \brief Implementation of Evans model for an elastic cytoskeleton & bilayer.

*/
#include <iostream>
#include "EvansElastic.h"

namespace voom 
{
  void EvansElastic::updateState(bool f0, bool f1, bool f2)
  {
    // first call updateState from the SCElastic to get bending
    // contributions
    SCElastic::updateState(f0, f1, f2);

    // OK, now add the stretching response
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }

    Tensor3D C(0.0);
    C = trans(F)*F;

    _trC = C(0,0)+C(1,1)+C(2,2);
    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    _J = sqrt((_trC*_trC - trCSquare)/2.0); 


    Tensor3D P(0.0);
    P = (_mu/_J)*F 
      + (1.0/_J)*(_kS*(_J-1) - 0.5*_mu*_trC/sqr(_J))*(_trC*F - F*C);
    
    if( f0 ) {
      // compute strain energy
      _Ws = (_kS*sqr(_J-1.0)/2.0 + 0.5*_mu*(_trC/_J-2.0))/_J;
      _W += _Ws;
    }
    
    if( f1 ) {
      //
      // -------------- stress resultants --------------
      for(int alpha=0; alpha<2; alpha++) 
	_n(alpha) += P*refDual(alpha)/_J;

      _n(2) = 0.0, 0.0, 0.0;
    }
    
    return;
  }

  const Tensor3D EvansElastic::DefGradient(){
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
   return F;
  }
  
  const Tensor3D & EvansElastic::cauchyStress() {
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
    Tensor3D C(0.0);
    C = trans(F)*F;

    _trC = C(0,0)+C(1,1)+C(2,2);
    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    _J = sqrt((_trC*_trC - trCSquare)/2.0); 
    Tensor3D P(0.0);
    P = (_mu/_J)*F 
      + (1.0/_J)*(_kS*(_J-1) - 0.5*_mu*_trC/sqr(_J))*(_trC*F - F*C);
   
    _cauchy=1/_J*P*trans(F);
    return _cauchy;
  }

  const std::vector<double > EvansElastic::invariants()
  {
    std::vector<double > invariants(2, 0.0);
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
    Tensor3D C(0.0);
    C = trans(F)*F;

    _trC = C(0,0)+C(1,1)+C(2,2);

    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
    _J = sqrt((_trC*_trC - trCSquare)/2.0);
     
    invariants[0] = _trC/_J;
    invariants[1] = _J;
    
    return invariants;
  }
}
