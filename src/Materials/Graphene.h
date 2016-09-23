// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Graphene_h__)
#define __Graphene_h__

#include<vector>
#include "SCElastic.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Class implementing the Foppl von Karman thin shell
    elasticity model.
  */
  class Graphene : public ShellMaterial<void*> {
  protected:
    
    // Covariant components of Green Strain
    Tensor2D _strain;

    // Covariant components of artificial reference Strain used for viscosity
    Tensor2D _Vstrain;

    // Contravariant components of 2nd P-K Stress
    Tensor2D _stress;

    //! Angle of the zig-zag direction relative to the x-axis
    double _theta;

  public:

    Graphene(double theta=0.0) 
    {
      _W = 0.0;
      Vector3D zero(0.0);
      _n(0) = zero;
      _n(1) = zero;
      _n(2) = zero;
      _m(0) = zero;
      _m(1) = zero;
      _stretch = 0.0;

      _theta = fmod(theta , M_PI/6.0);
      
      return;

    };

    void updateState(bool f0, bool f1, bool f2 );

    virtual double bendingEnergy() const {return 0.0;}
    virtual double stretchingEnergy() const {return _W;}

    const Tensor2D & strain() const {return _strain;}
    const Tensor2D & stress() const {return _stress;}

    void viscousStep() {
      _Vstrain = _strain;
    }
  };
  
} //namespace voom

#endif //  !defined(__Graphene_h__)
