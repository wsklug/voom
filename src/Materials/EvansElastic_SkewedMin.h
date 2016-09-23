// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                         Ankush Aggarwal
//                          Luigi Perotti
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file EvansElastic_Skewedmin.h

  \brief Material class based on EvansElastic_Skewed with the difference that the shear value is a dof

*/

#if !defined(__EvansElastic_SkewedMin_h__)
#define __EvansElastic_SkewedMin_h__

#include<blitz/array.h>
 
 
#include<vector>
#include "voom.h"
#include "../Node/Node.h"
#include "SCElastic.h"
#include "VoomMath.h"

using namespace std;

namespace voom
{

  /*!  Class implementing the Evans model for an elastic lipid bilayer & cytoskeleton. */
  
  class EvansElastic_SkewedMin : public ShellMaterial<void*> {

  protected:
    //! shear modulus
    double _mu; 

    //! stretching modulus
    double _kS;

    //! stretching energy
    double _Ws;

    //! second equilibrium value of shear and it's direction
    double _theta;
    double _eta;
    double _shearForce;      // Force conjugated to the shear magnitude
    double _directionForce;  // Force conjugated to the shear direction

  public:
    EvansElastic_SkewedMin(double mu, const double kS, const double eta, const double theta) 
      : _mu(mu), _kS(kS), _Ws(0.0), _eta(eta), _theta(theta), _shearForce(0.0), _directionForce(0.0) {};

    virtual void updateState(bool f0, bool f1, bool f2 );

    double shearModulus() const {return _mu;}
    void setShearModulus(double new_mu) { _mu = new_mu;}
    double stretchingModulus() const {return _kS;}
    void setStretchModulus(double new_kS) { _kS = new_kS;}
    void setEta(double new_eta) { _eta = new_eta;}
    void setTheta(double new_theta) { _theta = new_theta;}
    double shearForce() const {return _shearForce;}
    double directionForce() const {return _directionForce;}
    double stretchingEnergy() const {return _Ws;}
    double bendingEnergy() const {return 0.0;}
    const Tensor3D DefGradient();
    const Tensor3D & cauchyStress();
    const std::vector<double > invariants();

  };
  
} //namespace voom

#endif //  !defined(__EvansElastic_SkewedMin_h__)
