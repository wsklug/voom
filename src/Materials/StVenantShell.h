// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2013 All Rights Reserved
//
//----------------------------------------------------------------------
//

/*! 
  \file StVenantShell.h

  \brief Interface for St. Venant-Kirchhoff model for
  an elastic thin shell.

*/

#if !defined(__StVenantShell_h__)
#define __StVenantShell_h__

#include "ShellMaterial.h"
#include "VoomMath.h"


namespace voom
{

  /*!  Class implementing the simple Helfrich spontaneous curvature
    model for an elastic lipid bilayer.
  */
  
  class StVenantShell : public ShellMaterial<void*>
  {
  protected:

    double _E;
    double _D;
    double _nu;
    double _h;
    double _lambda;
    double _mu;
		
    Tensor2D _strain;    // strain tensor
    Tensor2D _stress;    // stress tensor
    Tensor2D _curvature; // curvature strain tensor
    Tensor2D _moments;   // moment resultant tensor

    double _H; // mean curvature
    double _K; // gaussian curvature

    double _Wb;
    double _Ws;

  public:

    StVenantShell(const double E, const double nu, const double h)
      : _E(E), _nu(nu), _h(h), _curvature(0.0) {
      _D = E*h/(12.0*(1-sqr(nu)));
      _lambda = nu*E*h/(1-sqr(nu));
      _mu = E*h/(2*(1+nu));
      
      _H = _K = _Wb = _Ws = _W;

    }

    virtual void updateState(bool f0, bool f1, bool f2 );

    double meanCurvature() const {return _H;}
    double gaussianCurvature() const {return _K;}

    double YoungsModulus() const { return _E; }
    double PoissonRatio() const { return _nu; }
    double bendingModulus() const { return _D; }
    double thickness() const { return _h; }

    virtual double bendingEnergy() const {return _Wb;}
    virtual double stretchingEnergy() const {return _Ws;}


  };
  
} //namespace voom

#endif //  !defined(__StVenantShell_h__)
