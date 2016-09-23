// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file ModEvansElastic.h

  \brief Material class for modeling conformational change using double welled potential.

*/

#if !defined(__ModEvansElastic_h__)
#define __ModEvansElastic_h__

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
  
  class ModEvansElastic : public SCElastic {

  protected:
    //! shear modulus
    double _mu; 

    //! stretching modulus
    double _kS;

    //! first invariant of Rt. Cauchy-Green deformation tensor
    double _trC, _trChat;

    //! second invariant of Rt. Cauchy-Green deformation tensor
    double _J, _Jhat;

    //! stretching energy
    double _Ws;

    //! second equilibrium value of shear and it's direction
    double _eta, _theta;

    //! material constants added to the energy terms
    double _c1, _c2;

    //! inverse of the deformation gradient of sheared equilibrium state
    Tensor3D _invhatF;

  public:
//need to add shear direction in the argument list
    ModEvansElastic(const double kC, const double kG, const double C0, double mu, 
		 const double kS, const double eta, const double theta, const double c1, const double c2 ) 
      : SCElastic(kC, kG, C0), _mu(mu), _kS(kS), _J(0.0), _trC(0.0), _Ws(0.0), _eta(eta), _theta(theta), _c1(c1), _c2(c2)
    {//see if it makes sense to add the formula for invhatF and do that
      _invhatF(0,2)=0.; _invhatF(1,2)=0.; _invhatF(2,0)=0.; _invhatF(2,1)=0.;
      _invhatF(0,0)=1+_eta*sin(_theta)*cos(_theta);
      _invhatF(0,1)=-_eta*cos(_theta)*cos(_theta);
      _invhatF(1,0)=_eta*sin(_theta)*sin(_theta);
      _invhatF(1,1)=1-_eta*sin(_theta)*cos(_theta);
      _invhatF(2,2)=0;
};

    virtual void updateState(bool f0, bool f1, bool f2 );
    //check consistency
    //void ConsistencyTest();

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

#endif //  !defined(__ModEvansElastic_h__)
