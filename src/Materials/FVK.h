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
    
  public:

    FVK(const double kC, const double kG, const double C0, const double E,
	const double nu) : SCElastic(kC, kG, C0), _E(E), _nu(nu)
    {
      _lambda = E*nu/( (1.0-nu)*(1.0+nu) );
      _mu = 0.5*E/(1.0+nu);
    };

    void updateState(bool f0, bool f1, bool f2 );

    double youngsModulus() const { return _E; }

    double poissonRatio() const { return _nu; }

    void setSpontaneousCurvature( double C0 ) { _C0 = C0; }

    virtual double bendingEnergy() const {return _W - _Ws;}
    virtual double stretchingEnergy() const {return _Ws;}
  };
  
} //namespace voom

#endif //  !defined(__FVK_h__)
