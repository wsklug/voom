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
  \file EvansElastic_Stretch.h

  \brief Material class based on EvansElastic_Skewed. 
  A pre-stretch along the principal directions is imposed instead of shear

*/

#if !defined(__EvansElastic_Stretch_h__)
#define __EvansElastic_Stretch_h__

#include "voom.h"
#include<blitz/array.h>
#include<vector>
#include "../Node/Node.h"
#include "SCElastic.h"
#include "VoomMath.h"

using namespace std;

namespace voom
{

  /*!  Class implementing the Evans model for an elastic lipid bilayer & cytoskeleton. */
  
  class EvansElastic_Stretch : public ShellMaterial<void*> {

  protected:
    //! shear modulus
    double _mu; 

    //! stretching modulus
    double _kS;

    //! stretching energy
    double _Ws;

    //! conformational energy
    double _Wc;

    //! second equilibrium value of stretch and it's direction
    double _theta;
    double _eta;
    double _phi;
    double _stretchForce;    // Force conjugated to the stretch magnitude
    double _directionForce;  // Force conjugated to the stretch direction
    double _phiForce;        // Force conjugated to the soft mode direction
    int _WcType;             // Type of conformational energy
    double* _WcConst;        // Multiplicative factor for conformational energy
    double c1, c2, c3, c4;   // Constant for W6 (to be computed only once based on _WcConst

  public:
    EvansElastic_Stretch(double mu, const double kS, const double eta, const double theta,  const int WcType = 0, double* WcConst = NULL) 
      : _mu(mu), _kS(kS), _Ws(0.0), _Wc(0.0), _eta(eta), _theta(theta), _phi(0.0), _stretchForce(0.0),  _directionForce(0.0), _phiForce(0.0), _WcType(WcType), _WcConst(WcConst) {

      if (_WcType == 6) {
	double eta0 = _WcConst[3];
	double eps = _WcConst[1];
	c1 = 0.5*_WcConst[2]*(pow(pow(eps, 2.0) - 1.0, 2.0)/(2.0*pow(eps, 3.0)*eta0) + ((2.0*eps - 1.0)*pow(eps + 1.0, 2.0))/(4.0*pow(eps, 3.0)*pow(eta0, 2.0)) - (eta0*pow(eps - 1.0, 2.0)*(4.0*eps - eta0 - 2.0*eps*eta0 + 2.0*pow(eps, 2.0) + 2.0))/(4.0*pow(eps, 3.0)));
	c2 = 0.5*_WcConst[2]*((3.0*pow(eta0, 4.0) - 2.0*pow(eta0, 3.0) + 2.0*eta0 - 3.0)/(4.0*eps*pow(eta0, 2.0)) - (pow(eta0, 2.0) + 1.0)/eta0 - (3.0*pow(eta0 - 1.0, 3.0)*(eta0 + 1.0))/(4.0*pow(eps, 3.0)*pow(eta0, 2.0)));
	c3 = 0.5*_WcConst[2]*((3.0*pow(eta0 - 1.0, 3.0)*(eta0 + 1.0))/(4.0*pow(eps, 3.0)*pow(eta0, 2.0)) - (pow(eta0, 2.0) - 1)/(2.0*eps*eta0) + 1.0);
	c4 = -0.5*_WcConst[2]*((pow(eta0 - 1.0, 3.0)*(eta0 + 1.0))/(4.0*pow(eps, 3.0)*pow(eta0, 2.0)));
      }
      else {
	c1 = 0.0;
	c2 = 0.0;
	c3 = 0.0;
	c4 = 0.0;
      }
};

    virtual void updateState(bool f0, bool f1, bool f2 );

    double shearModulus() const {return _mu;}
    void setShearModulus(double new_mu) { _mu = new_mu;}
    double stretchingModulus() const {return _kS;}
    void setStretchModulus(double new_kS) { _kS = new_kS;}
    void setEta(double new_eta) { _eta = new_eta;}
    void setTheta(double new_theta) { _theta = new_theta;}
    void setPhi(double new_phi) { _phi = new_phi;}
    double stretchForce() const {return _stretchForce;}
    double directionForce() const {return _directionForce;}
    double phiForce() const {return _phiForce;}
    double stretchingEnergy() const {return _Ws;}
    double conformationalEnergy() const {return _Wc;}
    double bendingEnergy() const {return 0.0;}
    const Tensor3D DefGradient();
    const Tensor3D & cauchyStress();
    const std::vector<double > invariants();

  };
  
} //namespace voom

#endif //  !defined(__EvansElastic_Stretch_h__)
