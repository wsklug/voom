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
    FVK::updateState(f0||f1, f1, f2); //claculate _n & _m, _W 

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
    
    double g=0.0, dg=0.0;
    if(_formulation==0) {
      // g is a simple quartic with curvature g0 at eta=0 and g(1)=g(0)+_deltag
      double eta2=_eta*_eta;
      double eta3=eta2*_eta;
      double eta4=eta3*_eta;
      g  = _g0*(eta4-2.0*eta2) - 2.0*(_deltag+_g0)*(eta3-1.5*eta2);
      dg = (_g0*(4.0*_eta - 2.0) - 6.0*_deltag)*(eta2-_eta); 
    } else {
      // g is a two piece function which has different curvature at two minima
      double g1=_g0-32.0*_deltag;

      if (_eta < 0.5){
	g  = _g0*_eta*_eta*(_eta*_eta - 2.0*_eta + 1.0)/2.0;
	dg = _g0*_eta*(2.0*_eta*_eta - 3.0*_eta + 1.0);
      }
      
      else {
	g  = g1*_eta*_eta*(_eta*_eta - 2.0*_eta + 1.0)/2.0 + _deltag;
	dg = g1*_eta*(2.0*_eta*_eta - 3.0*_eta + 1.0);
      }
    }

    // subtract a constant chemical potential term
    g  -= _muA*_eta;
    dg -= _muA;

    double strainEnergy=0.0;
    if(f1) strainEnergy = _W;

    if (f0){
      _W *= _eta;

      _Ws *= _eta;

      _W += g*jacobian;

      //_W += -_gammaZero*_eta*traceStrain*jacobian;
      
      for (int alpha = 0; alpha < 2; alpha++){
	for (int beta =0; beta <2; beta++){

	  _W += 0.5*_Gamma*_deta(alpha)*_deta(beta)*refMetricInv(alpha, beta)*jacobian;

	}
      }

    }



    if (f1){
      for (int alpha = 0; alpha < 2; alpha++){
// 	_n(alpha) = 0.0, 0.0, 0.0;
	_n(alpha) *= _eta;

// 	for (int mu = 0; mu < 2; mu++){
	  //_n(alpha) += -_gammaZero*_eta*refMetricInv(alpha, mu)*basis(mu)*jacobian;

	// This only when integrating over deformed shape
// 	  for (int beta = 0; beta < 2; beta ++){
// 	    _n(alpha) += - _Gamma*refMetricInv(alpha, mu)*_deta(mu)*_deta(beta)*dual(beta)*jacobian;

// 	  }
// 	}
      }

      _n(2)= 0.0, 0.0, 0.0;


//       _mu = (dg - _gammaZero*traceStrain)*jacobian;
      _mu = dg*jacobian;

      _mu += strainEnergy;

      for (int beta = 0; beta < 2; beta++){
	_lambda(beta) = 0.0;

	for (int alpha = 0; alpha < 2; alpha++){
	  _lambda(beta) += _Gamma*refMetricInv(alpha, beta)*_deta(alpha)*jacobian;

	}
      }
    }
    return;
  }
}
