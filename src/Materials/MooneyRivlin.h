// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  Ankush Aggarwal, William S. Klug
//                University of California Los Angeles
//                    (C) 2011 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file MooneyRivlin.h

  \brief Interface for a compressible MooneyRivlin finite deformation
  hyperelasticity model.

*/

#ifndef _MOONEYRIVLIN_H_
#define _MOONEYRIVLIN_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {
  
  class MooneyRivlin : public Material
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    MooneyRivlin() { _init(0.0,0.0,0.0); }
    //! Construct material from density, Young's modulus, and Poisson's ratio
    MooneyRivlin(double C10, double C01, double D1) {
      _init(C10, C01, D1);
    }
    //! Destructor
    virtual ~MooneyRivlin() {}
    //! Copy constructor
    MooneyRivlin(const MooneyRivlin &);
    
    // Accessors/mutators:
    //! Returns density
   // inline double massDensity();
    //! Calculates longitudinal wave speed
   // inline double longitudinalWaveSpeed();
    
    // Operators
    
    // General methods:
    //! Based on new deformation gradient tensor, F, calculates state of material (strain energy density, first Piola-Kirchhoff stress tensor)
    void updateState(bool f0, bool f1, bool f2);
    double vonMisesStress() const;
    // Tests:
    //! Consistency test
    void ConsistencyTest();
   // static void MFITest();
   // static void IsotropyTest();
    void setModulus(double C10, double C01, double D1){_C10=C10; _C01=C01; _D1=D1;};
    //double getYoungsModulus(){return _E;};

    //copy function is pure virtual
    Material * copy() {};
    
  private:
    
    // Members:
    //! Moduli
    double _C10,_C01,_D1;
    
    //! Initializes material properties and zeros out deformation/stress tensors
    void _init(double rho, double E, double nu);
  };
  
}
#endif // _MOONEYRIVLIN_H_
