// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file CompNeoHookean.h

  \brief Interface for a compressible Neo-Hookean finite deformation
  hyperelasticity model.

*/

#ifndef _COMPNEOHOOKEAN_H_
#define _COMPNEOHOOKEAN_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {
  
  class CompNeoHookean : public Material
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    CompNeoHookean() { _init(0.0,0.0,0.0); }
    //! Construct material from density, Young's modulus, and Poisson's ratio
    CompNeoHookean(double rho, double E, double nu) {
      _init(rho, E, nu);
    }
    //! Destructor
    virtual ~CompNeoHookean() {}
    //! Copy constructor
    CompNeoHookean(const CompNeoHookean &);
    
    // Accessors/mutators:
    //! Returns density
    inline double massDensity();
    //! Calculates longitudinal wave speed
    inline double longitudinalWaveSpeed();
    
    // Operators
    
    // General methods:
    //! Based on new deformation gradient tensor, F, calculates state of material (strain energy density, first Piola-Kirchhoff stress tensor)
    void updateState(bool f0, bool f1, bool f2);
    double vonMisesStress() const;
    // Tests:
    //! Consistency test
    void ConsistencyTest();
    static void MFITest();
    static void IsotropyTest();
    void setYoungsModulus(double E){_E=E;};
    double getYoungsModulus(){return _E;};

    //copy function is pure virtual
    Material * copy() {};
    
  private:
    
    // Members:
    //! Density
    double _rho;
    //! Young's Modulus
    double _E;
    //! Poisson's ratio
    double _nu;
    
    //! Initializes material properties and zeros out deformation/stress tensors
    void _init(double rho, double E, double nu);
  };
  
}
#endif // _COMPNEOHOOKEAN_H_
