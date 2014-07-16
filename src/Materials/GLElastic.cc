// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log: EvansElastic.cc,v $
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
  \file GLElastic.h

  \brief Interface for structural phase transformation based on Ginzburg-Landau theory.

*/
#include <iostream>
#include "GLElastic.h"

namespace voom 
{
  void GLElastic::updateState(bool f0, bool f1, bool f2)
  {
    FVK::updateState(f0, f1, f2); //claculate _n & _m, _W 

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();
    const BasisVectors & dPartials = _deformedGeometry.dPartials();
    
    const BasisVectors & rbasis = _referenceGeometry.a();
    const BasisVectors & rdual = _referenceGeometry.aDual();
    
    _H = - 0.5 * ( dot(dual(0), dPartials(0)) + dot(dual(1), dPartials(1)) );
    
    const Tensor2D& metricInv = 
      _deformedGeometry.metricTensorInverse();
    const Tensor2D& refMetricInv =
      _referenceGeometry.metricTensorInverse();

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();

    double jacobian = _referenceGeometry.metric()/_deformedGeometry.metric();

    Tensor2D strain(0.0);

    strain = 0.5*(metric-refMetric);
//     if (strain(0,1) != strain(1,0)){
//       double xx=strain(0,1);
//       std::cout<<"xx="<<xx<<std::endl;
//       double yy=strain(1,0);
//       std::cout<<"yy="<<yy<<std::endl;
//       std::cout<<"xx-yy="<<xx-yy<<std::endl;
//     }
    assert(std::abs(strain(0,1)-strain(1,0))<1.0e-8);
    if (std::abs(strain(0,1) - strain(1,0))>1.0e-8){
      double xx=strain(0,1);
      std::cout<<"xx="<<xx<<std::endl;
      double yy=strain(1,0);
      std::cout<<"yy="<<yy<<std::endl;
      std::cout<<"xx-yy="<<xx-yy<<std::endl;
    }
    double traceStrain = 0.0;
    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	traceStrain   += refMetricInv(alpha,beta)*strain(alpha,beta);
      }
    }

		
    double twoHminusC0 = 2.0*_H - _C0;


//     double g = _gDoublePrime*_eta*_eta*(_eta*_eta - 2.0*_eta + _c)/2.0/_c;
//     double dg = _gDoublePrime*_eta*(2.0*_eta*_eta - 3.0*_eta + _c)/_c;
    
    double g, dg;
    // g is now a two piece function which has different curvature at two minima
    if (_eta < 0.5){
      g  = _gDoublePrime*_eta*_eta*(_eta*_eta - 2.0*_eta + 1.0)/2.0;
      dg = _gDoublePrime*_eta*(2.0*_eta*_eta - 3.0*_eta + 1.0);
    }

    else {
      g  = _gDoublePrime*_eta*_eta*(_eta*_eta - 2.0*_eta + 1.0)/2.0/_factor + _gDoublePrime*(1.0-1.0/_factor)/32.0;
      dg = _gDoublePrime*_eta*(2.0*_eta*_eta - 3.0*_eta + 1.0)/_factor;
    }

    if (f0){
      _W += g*jacobian;

      _W += -_gammaZero*_eta*traceStrain*jacobian;
      
      for (int alpha = 0; alpha < 2; alpha++){
	for (int beta =0; beta <2; beta++){

	  _W += _Gamma/2.0*_dEta(alpha)*_dEta(beta)*metricInv(alpha, beta)*jacobian;

	}
      }

    }



    if (f1){
      for (int alpha = 0; alpha < 2; alpha++){
	_nEta(alpha) = 0.0, 0.0, 0.0;

	for (int mu = 0; mu < 2; mu++){
	  _nEta(alpha) += -_gammaZero*_eta*refMetricInv(alpha, mu)*basis(mu)*jacobian;

	  for (int beta = 0; beta < 2; beta ++){
	    _nEta(alpha) += - _Gamma*metricInv(alpha, mu)*_dEta(mu)*_dEta(beta)*dual(beta)*jacobian;

	  }
	}
      }

      _nEta(2)= 0.0, 0.0, 0.0;


      _muEta = (dg - _gammaZero*traceStrain)*jacobian;

      for (int beta = 0; beta < 2; beta++){
	_muDEta(beta) = 0.0;

	for (int alpha = 0; alpha < 2; alpha++){
	  _muDEta(beta) += _Gamma*metricInv(alpha, beta)*_dEta(alpha)*jacobian;

	}
      }
    }
    return;
  }
}
