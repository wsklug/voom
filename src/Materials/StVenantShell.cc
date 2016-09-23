// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2013 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file StVenantShell.cc

  \brief Implementation of the simple Helfrich spontaneous curvature
  model for an elastic lipid bilayer.

*/
#include <iostream>
#include "StVenantShell.h"

namespace voom 
{
  void StVenantShell::updateState(bool f0, bool f1, bool f2)
  {

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;

    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();
    const BasisVectors & dPartials = _deformedGeometry.dPartials();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    const BasisVectors & refDPartials = _referenceGeometry.dPartials();

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();
    const Tensor2D& refMetricInv = _referenceGeometry.metricTensorInverse();

    //
    // Compute middle surface strains
    //
    _strain = 0.5*(metric-refMetric);

    //
    // Compute curvature difference tensor
    //
    _curvature = Tensor2D(0.0);

    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	_curvature(alpha,beta) -= dot( basis(alpha), dPartials(beta) )
	                        - dot( refBasis(alpha), refDPartials(beta) );
      }
    }

    Tensor2D strainDual(0.0);
    Tensor2D curvatureDual(0.0);
    strainDual = refMetricInv*_strain*tvmet::trans(refMetricInv);
    curvatureDual = refMetricInv*_curvature*tvmet::trans(refMetricInv);

    double traceStrain = 0.0;
    double strainSquared = 0.0;
    double traceCurvature = 0.0;
    double curvatureSquared = 0.0;
    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	traceStrain   += refMetricInv(alpha,beta)*_strain(alpha,beta);
	strainSquared += strainDual(alpha,beta)*_strain(alpha,beta);
	traceCurvature   += refMetricInv(alpha,beta)*_curvature(alpha,beta);
	curvatureSquared += curvatureDual(alpha,beta)*_curvature(alpha,beta);
      }
    }

    // contravariant components of 2nd P-K stress
    _stress = _lambda*traceStrain*refMetricInv + 2.0*_mu*strainDual;

    // bending moments

    _moments = _D*(  _nu*traceCurvature*refMetricInv 
		     + 2.0*(1.0-_nu)*curvatureDual);

    // one-half Trace of curvature difference tensor
    _H = 0.5*traceCurvature;

    _K =  (  _curvature(0,0) * _curvature(1,1) 
	   - _curvature(0,1) * _curvature(1,0) )/_referenceGeometry.metric();

    double jacobian = _referenceGeometry.metric()/_deformedGeometry.metric();

    // compute strain energy

    if( f0 ) {

      _Ws = (0.5*_lambda*traceStrain*traceStrain + _mu*strainSquared)*jacobian;

      _Wb = 0.5*_D*( _nu*sqr(2.0*_H) + (1.0-_nu)*curvatureSquared );
      
      _W = _Ws + _Wb;
    }
		
    if( f1 ) {
      //
      // -------------- stress resultants --------------
      for(int alpha=0; alpha<2; alpha++) {
	_n(alpha) = 0.0,0.0,0.0;
	for(int beta=0; beta<2; beta++) {
	  _n(alpha) += _stress(alpha,beta)*basis(beta)*jacobian;

	  _n(alpha) -= _moments(alpha,beta)*dPartials(beta)*jacobian;
	}
      }

      // 
      _n(2) = 0.0, 0.0, 0.0;
      //
      // ------------------- moment resultants ------------------
      for(int alpha=0; alpha<2; alpha++) {
	_m(alpha) = 0.0,0.0,0.0;
	for(int beta=0; beta<2; beta++) {
	  _m(alpha) -= _moments(alpha,beta)*basis(beta)*jacobian;
	}
      }
      
    }

    return;
  }
}
