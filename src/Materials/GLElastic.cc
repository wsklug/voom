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
    
    const Tensor2D& metricInv = 
      _deformedGeometry.metricTensorInverse();
    const Tensor2D& refMetricInv =
      _referenceGeometry.metricTensorInverse();

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();

    // Previously we'd considered the GL energy as defined per unit
    // reference area.  Since shell materials are defined per unit
    // deformed area (who knows why), we needed a Jacobian, the ratio
    // between referenced and deformed areas.  Now we have changed
    // this energy to be per unit deformed area, so we no longer need
    // the jacobian.
    
    //double jacobian = _referenceGeometry.metric()/_deformedGeometry.metric();
		
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

    double W_GL = 0.0;
    if(f0||f1) {
      // now add GL energy
      W_GL += g;//*jacobian;

      for (int alpha = 0; alpha < 2; alpha++){
	for (int beta =0; beta <2; beta++){

	  // if gradient is defined on reference configuration
	  //W_GL += 0.5*_Gamma*_deta(alpha)*_deta(beta)*refMetricInv(alpha, beta);

	  // if gradient is defined on deformed configuration
	  W_GL += 0.5*_Gamma*_deta(alpha)*_deta(beta)*metricInv(alpha, beta);
	}
      }
    }

    if (f0){
      // scale elastic terms by order parameter for linear coupling
      _W *= _eta;

      _Ws *= _eta;

      _stress *= _eta;

      _W += W_GL;
    }


    // Compute stress & chemical potential resultants

    if (f1){
      // stress resultants
      for (int alpha = 0; alpha < 2; alpha++){
	// scale elastic term by order parameter for linear coupling
	_n(alpha) *= _eta;

	// variation of deformed area metric
	_n(alpha) += W_GL*dual(alpha);
      }

      _n(2)= 0.0, 0.0, 0.0;


      // chemical potential resultants
      _mu = strainEnergy;
      _mu += dg;


      // chemical potential gradient resultants
      for (int alpha = 0; alpha < 2; alpha++){
	_lambda(alpha) = 0.0;

	for (int beta = 0; beta < 2; beta++){
	  // if gradient is defined on reference configuration
	  //_lambda(alpha) += _Gamma*refMetricInv(alpha, beta)*_deta(beta);
	  
	  // if gradient is defined on deformed configuration
	  _lambda(alpha) += _Gamma*metricInv(alpha, beta)*_deta(beta);

	  for (int mu=0; mu<2; mu++){
	    _n(alpha) += 
	      - _Gamma*metricInv(alpha, mu)*_deta(mu)*_deta(beta)*dual(beta);
	  }

	}
      }
    }
    return;
  }
}
