// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file EvansElastic.h

  \brief Interface for Evans model for an elastic cytoskeleton & bilayer.

*/

#if !defined(__EvansElastic_h__)
#define __EvansElastic_h__

#include<blitz/array.h>
 
 
#include<vector>
#include "voom.h"
#include "SCElastic.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Class implementing the Evans 
    model for an elastic lipid bilayer & cytoskeleton.
  */
  
  class EvansElastic : public SCElastic {

  protected:
    //! shear modulus
    double _mu; 

    //! stretching modulus
    double _kS;

    //! first invariant of Rt. Cauchy-Green deformation tensor
    double _trC;

    //! second invariant of Rt. Cauchy-Green deformation tensor
    double _J;

    //! stretching energy
    double _Ws;

  public:

    EvansElastic(const double kC, const double kG, const double C0, double mu, 
		 const double kS ) 
      : SCElastic(kC, kG, C0), _mu(mu), _kS(kS), _J(0.0), _trC(0.0), _Ws(0.0)
    {};

    virtual void updateState(bool f0, bool f1, bool f2 );

    double shearModulus() const {return _mu;}
    double stretchingModulus() const {return _kS;}
    double J() const {return _J;}
    double trC() const {return _trC;}

    void setShearModulus(double new_mu) {
      _mu = new_mu;
    }

    void setStretchModulus(double new_kS) {
      _kS = new_kS;
    }

    virtual double bendingEnergy() const {return _W - _Ws;}
    virtual double stretchingEnergy() const {return _Ws;}
  };
  
} //namespace voom

#endif //  !defined(__EvansElastic_h__)
