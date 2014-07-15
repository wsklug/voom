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
// Revision 1.4  2005/05/23 17:36:54  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file Hookean2D.cc

  \brief Concrete class for an isotropic, 2-D (plane stress/strain),
  small-strain, Hookean material object.

*/

#include "Hookean2D.h"

namespace voom
{

/*! Concrete class for an isotropic, 2-D (plane stress/strain),
  small-strain, Hookean material object.
*/
  Hookean2D::Hookean2D( Option o, double E, double nu ) {
    _option = o;
    _E = E;
    _nu = nu;
    _lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    _mu = 0.5*E/(1.0+nu);
    
    // initialize stuff
    _w = 0.0;
    _epsilon = 0.0, 0.0, 0.0;
    _sigma = 0.0, 0.0, 0.0;

    // matrix of elastic moduli is constant, so you may want to
    // compute it now in the constructor so that you never have to
    // again.
    double c11=0.0,c12=0.0;
    if( _option == planeStrain ) {
      c11 = _lambda + 2.0*_mu;
      c12 = _lambda;
    } else {
      c11 = _E/(1.0-_nu*_nu);
      c12 = _E*_nu/(1.0-_nu*_nu);
    }
    _c = c11, c12, 0.0,
         c12, c11, 0.0,
         0.0, 0.0, _mu;
  } 

  void Hookean2D::updateState(bool f0, bool f1, bool f2) {

    if(f0 || f1) {
      // Compute stress
      _sigma = _c*_epsilon;
    }
    if(f0) {
      // Compute strain energy
      _w = 0.5* dot(_sigma,_epsilon); 
    }

    return;
  }
} //namespace voom

