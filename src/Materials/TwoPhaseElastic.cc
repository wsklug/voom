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
  \file TwoPhaseElastic.cc

  \brief Implementation for two phase Evans model for an elastic cytoskeleton & bilayer.

*/
#include <iostream>
#include "TwoPhaseElastic.h"

namespace voom 
{
  void TwoPhaseElastic::updateState(bool f0, bool f1, bool f2)
  {

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();
    const BasisVectors & dPartials = _deformedGeometry.dPartials();
    
    const BasisVectors & rbasis = _referenceGeometry.a();
    const BasisVectors & rdual = _referenceGeometry.aDual();
    
    _H = - 0.5 * ( dot(dual(0), dPartials(0)) + dot(dual(1), dPartials(1)) );
    
    const Tensor2D& aInv = 
      _deformedGeometry.metricTensorInverse();

    Tensor2D mixedCurvatureTensor;
    // b^alpha_beta = mixedCurvatureTensor(alpha,beta)
    // first index is upper, second is lower
    mixedCurvatureTensor(0,0) = -dot(dual(0), dPartials(0));
    mixedCurvatureTensor(0,1) = -dot(dual(0), dPartials(1));
    mixedCurvatureTensor(1,0) = -dot(dual(1), dPartials(0));
    mixedCurvatureTensor(1,1) = -dot(dual(1), dPartials(1));

    _K =  (  mixedCurvatureTensor(0,0) * mixedCurvatureTensor(1,1) 
	     - mixedCurvatureTensor(0,1) * mixedCurvatureTensor(1,0) );
		
    double twoHminusC0 = 2.0*_H - _C0;

    setkCkGC0();    

    EvansElastic::updateState(f0, f1, f2); //claculate _n & _m, _W & _WN

    if (f0){
      _W += _h*16.0*_C*_C*(_C-1.0)*(_C-1.0);
      
      for (int alpha = 0; alpha < 2; alpha++){
	for (int beta =0; beta <2; beta++){

	  _W += _epsilon*_dC(alpha)*_dC(beta)*aInv(alpha, beta);

	}
      }
      }


    if (f1){
      for (int alpha = 0; alpha < 2; alpha++){

	_nC(alpha) = _h*16.0*_C*_C*(_C-1.0)*(_C-1.0)*dual(alpha);

	for (int mu = 0; mu < 2; mu++){
	  for (int beta = 0; beta < 2; beta ++){

	    _nC(alpha) += _epsilon*aInv(mu, beta)*_dC(mu)*_dC(beta)*dual(alpha) - 2.0*_epsilon*aInv(alpha, mu)*_dC(mu)*_dC(beta)*dual(beta);

	  }
	}
      }

      _nC(2)= 0.0, 0.0, 0.0;


      _muC = -_kC*twoHminusC0*(_C01-_C02)+0.5*twoHminusC0*twoHminusC0*(_kC1-_kC2)+_K*(_kG1-_kG2)+_h*32.0*_C*(_C-1.0)*(2*_C-1.0);//l+h*psi'

      for (int beta = 0; beta < 2; beta++){
	_muDC(beta) = 0.0;

	for (int alpha = 0; alpha < 2; alpha++){
	  _muDC(beta) += 2.0*_epsilon*aInv(alpha, beta)*_dC(alpha);

	}
      }
    }
    return;
  }
}
