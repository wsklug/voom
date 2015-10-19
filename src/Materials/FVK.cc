// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/10/21 00:32:54  klug
// Initial checkin.
//
//
//----------------------------------------------------------------------

#include <iostream>
#include "FVK.h"

namespace voom 
{
  void FVK::updateState(bool f0, bool f1, bool f2)
  {
    SCElastic::updateState(f0, f1, f2);

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;

    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & refBasis = _referenceGeometry.a();

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();
    const Tensor2D& refMetricInv = _referenceGeometry.metricTensorInverse();

    Tensor2D strain;
    strain = 0.5*(metric-refMetric);

    Tensor2D strainDual(0.0);
    strainDual = refMetricInv*strain*tvmet::trans(refMetricInv);

    double traceStrain = 0.0;
    double strainSquared = 0.0;
    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	traceStrain   += refMetricInv(alpha,beta)*strain(alpha,beta);
	strainSquared += strainDual(alpha,beta)*strain(alpha,beta);
      }
    }

    double jacobian = _referenceGeometry.metric()/_deformedGeometry.metric();

    if( f0 ) {
      // compute strain energy
      _Ws = (0.5*_lambda*traceStrain*traceStrain + _mu*strainSquared)*jacobian;
      _W += _Ws;
    }
		
    if( f1 ) {

      // stress resultants
      for(int alpha=0; alpha<2; alpha++) {
	for(int beta=0; beta<2; beta++) {
	  _n(alpha) += 
	    ( 2.0*_mu*strainDual(alpha,beta)
	      + _lambda*traceStrain*refMetricInv(alpha,beta) ) 
	    * basis(beta) * jacobian;
	}
      }

      _n(2) = 0.0, 0.0, 0.0;
    }

    return;
  }
  
  //! Returns principal strains of Green-Lagrange Strain tensor
  double FVK::getMaxStrain() const{

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();
    const Tensor2D& refMetricInv = _referenceGeometry.metricTensorInverse();

    Tensor2D strain;
    strain = 0.5*(metric-refMetric);
    //To calculate the eigen values using same ways as we would for a
    //Cartesian coordinates matrix we need to transform the strain
    //tensor as follows
    strain = refMetricInv*strain;

    //Calculate eigen-values of strain.
    double T = strain(0,0) + strain(1,1);
    double D = strain(0,0)*strain(1,1) - strain(0,1)*strain(1,0);    
    double L1 = (0.5*T) + sqrt(0.25*T*T - D);
    double L2 = (0.5*T) - sqrt(0.25*T*T - D);
    
    //Return higher of two eigen values
    return (L1 > L2)? L1: L2;
  }
}
