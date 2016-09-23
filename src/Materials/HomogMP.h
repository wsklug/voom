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
  \file HomogMP.h

  \brief Interface for an experimental material HomogMP which minimizes at the most probable state in Homogenization technique.

*/

#ifndef _HOMOGMP_H_
#define _HOMOGMP_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {
  
  class HomogMP : public Material
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    HomogMP() { _init(1.0,1.0,1.0); }
    //! Construct material from density, Young's modulus, and Poisson's ratio
    HomogMP(double C10, double C01, double D1) {
      _init(C10, C01, D1);
    }
    //! Destructor
    virtual ~HomogMP() {}
    //! Copy constructor
    HomogMP(const HomogMP &);
    
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
    void setMPinv(double I1MP, double I2MP, double JMP){_I1MP=I1MP; _I2MP=I2MP; _JMP=JMP;};
    //double getYoungsModulus(){return _E;};

    //copy function is pure virtual
    Material * copy() {};
    
  private:
    
    // Members:
    //! Moduli
    double _C10,_C01,_D1,_I1MP,_I2MP,_JMP;
    
    //! Initializes material properties and zeros out deformation/stress tensors
    void _init(double C10, double C01, double D1);
  };
  
}
#endif // _HOMOGMP_H_
