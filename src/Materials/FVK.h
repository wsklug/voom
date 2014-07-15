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

#if !defined(__FVK_h__)
#define __FVK_h__

#include<vector>
#include "SCElastic.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Class implementing the Foppl von Karman thin shell
    elasticity model.
  */
  class FVK : public SCElastic {
  protected:
    double _E;
    double _nu;
    double _lambda;
    double _mu;
    double _Ws;
    
    // Covariant components of Green Strain
    Tensor2D _strain;

    // Contravariant components of 2nd P-K Stress
    Tensor2D _stress;
    
  public:

    FVK(const double kC, const double kG, const double C0, const double E,
	const double nu) : SCElastic(kC, kG, C0), _E(E), _nu(nu)
    {
      _lambda = E*nu/( (1.0-nu)*(1.0+nu) );
      _mu = 0.5*E/(1.0+nu);

      _W = 0.0;
      Vector3D zero(0.0);
      _n(0) = zero;
      _n(1) = zero;
      _n(2) = zero;
      _nTC(0) = zero;
      _nTC(1) = zero;
      _nTC(2) = zero;
      _m(0) = zero;
      _m(1) = zero;
      _stretch = 0.0;

    };

    void updateState(bool f0, bool f1, bool f2 );

    double youngsModulus() const { return _E; }

    double poissonRatio() const { return _nu; }

    virtual double bendingEnergy() const {return _W - _Ws;}
    virtual double stretchingEnergy() const {return _Ws;}

    const Tensor2D & strain() const {return _strain;}
    const Tensor2D & stress() const {return _stress;}

  };
  
} //namespace voom

#endif //  !defined(__FVK_h__)
