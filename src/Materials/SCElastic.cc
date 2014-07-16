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
// Revision 1.11  2005/10/10 00:55:22  klug
// Cleaned up comments.
//
// Revision 1.10  2005/08/20 01:28:31  klug
// acinclude.m4
//
// Revision 1.9  2005/06/18 17:05:02  klug
// Added Gaussian curvature term to energy.
//
// Revision 1.8  2005/05/23 17:38:39  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file SCElastic.cc

  \brief Implementation of the simple Helfrich spontaneous curvature
  model for an elastic lipid bilayer.

*/
#include <iostream>
#include "SCElastic.h"

namespace voom 
{
  void SCElastic::updateState(bool f0, bool f1, bool f2)
  {

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;

    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();
    const BasisVectors & dPartials = _deformedGeometry.dPartials();

    _H = - 0.5 * ( dot(dual(0), dPartials(0)) + dot(dual(1), dPartials(1)) );

    const Tensor2D& aInv = 
      _deformedGeometry.metricTensorInverse();
		
    // 		Tensor2D curvatureTensor;
    // 		curvatureTensor(0,0) = -dot(basis(0), dPartials(0));
    // 		curvatureTensor(0,1) = -dot(basis(0), dPartials(1));
    // 		curvatureTensor(1,0) = -dot(basis(1), dPartials(0));
    // 		curvatureTensor(1,1) = -dot(basis(1), dPartials(1));	       

    Tensor2D mixedCurvatureTensor;
    // b^alpha_beta = mixedCurvatureTensor(alpha,beta)
    // first index is upper, second is lower
    mixedCurvatureTensor(0,0) = -dot(dual(0), dPartials(0));
    mixedCurvatureTensor(0,1) = -dot(dual(0), dPartials(1));
    mixedCurvatureTensor(1,0) = -dot(dual(1), dPartials(0));
    mixedCurvatureTensor(1,1) = -dot(dual(1), dPartials(1));

    _K =  (  mixedCurvatureTensor(0,0) * mixedCurvatureTensor(1,1) 
	     - mixedCurvatureTensor(0,1) * mixedCurvatureTensor(1,0) );

    //
    //  edited by Feng on Dec. 7th, 2004
    //  Original code:
    //  double twoHminusC0 = 2.0*_H - _C0;
    //  modified code:
    //  Reason: we are using positive spontaneous curvatures (_C0), but a sphere has
    //          negative mean curvature (_H), -1/R, I believe that the operatore should be "+"
    double twoHminusC0 = 2.0*_H /* + */ - _C0;

    if( f0 ) {
      // compute strain energy
      _W =  0.5 * _kC * twoHminusC0 * twoHminusC0;
      _W += _kG * _K;
    }
		
    if( f1 ) {
      //
      // -------------- stress resultants --------------
      _n(0) = _kC * twoHminusC0 * 
	( aInv(0,0) * dPartials(0) + 
	  aInv(0,1) * dPartials(1) )
	+ 0.5 * _kC * twoHminusC0 * twoHminusC0 * dual(0) ;
      _n(0) += _kG*_K*dual(0);

      _n(1) = _kC * twoHminusC0 * 
	( aInv(1,0) * dPartials(0) + 
	  aInv(1,1) * dPartials(1) )
	+ 0.5 * _kC * twoHminusC0 * twoHminusC0 * dual(1) ;
      _n(1) += _kG*_K*dual(1); 

      for(int beta=0; beta<2; beta++) {
	for(int alpha=0; alpha<2; alpha++) {
	  _n(beta) += _kG*2.0*_H*aInv(beta,alpha)*dPartials(alpha);
	  for(int nu=0; nu<2; nu++) {
	    _n(beta) -= _kG*dot(dual(beta),dPartials(nu)) * 
	      mixedCurvatureTensor(nu,alpha)*dual(alpha);
	  }
	}
      }
      // 
      _n(2) = 0.0, 0.0, 0.0;
      //
      // ------------------- moment resultants ------------------
      _m(0) = - _kC * twoHminusC0 * dual(0);
      _m(1) = - _kC * twoHminusC0 * dual(1);

      for(int beta=0; beta<2; beta++) {
	_m(beta) -= _kG*2.0*_H*dual(beta);
	for(int alpha=0; alpha<2; alpha++) {
	  _m(beta) += _kG*mixedCurvatureTensor(beta,alpha)*dual(alpha);
	}
      }
    }

    return;
  }
}
