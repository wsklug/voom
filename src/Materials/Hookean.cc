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
// Revision 1.3  2005/05/23 17:36:54  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file Hookean.cc

  \brief Concrete class for an isotropic, 3-D, small-strain, Hookean
  material object.

*/

#include "Hookean.h"

namespace voom
{

/*! Concrete class for an isotropic, 3-D, small-strain, Hookean
  material object.
*/
  Hookean::Hookean( double E, double nu ) {
    _E = E;
    _nu = nu;
    _lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    _mu = 0.5*E/(1.0+nu);

    // initialize stuff
    _w = 0.0;
    VoigtVector zeroVec(0.0);
    VoigtMatrix zeroMat(0.0);
    _epsilon = zeroVec;
    _sigma = zeroVec;
    _c = zeroMat;

    // matrix of elastic moduli is constant, so compute it now and
    // never again
    double c11 = _lambda + 2.0*_mu;
    double c12 = _lambda;
    _c = c11, c12, c12, 0.0, 0.0, 0.0,
         c12, c11, c12, 0.0, 0.0, 0.0,
         c12, c12, c11, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, _mu, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, _mu, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, _mu;
  } 

  void Hookean::updateState(bool f0, bool f1, bool f2) {

    if(f0 || f1) {
      // Compute stress
      _sigma = _c*_epsilon;
    }
    if(f0) {
      // Compute strain energy
      _w = 0.5*tvmet::dot(_sigma,_epsilon); 
   }

    return;
  }
} //namespace voom

