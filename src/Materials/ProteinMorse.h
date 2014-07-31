// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file ProteninMorse.h
*/

#ifndef _PROTEINMORSE_H_
#define _PROTEINMORSE_H_

#include "ProteinPotential.h"

using namespace std;

namespace voom {
  
  class ProteinMorse: public ProteinPotential
  {
  public:
    // Constructors/destructors:
    //! Default constructor
    ProteinMorse(): _epsilon(0.0), _sigma(0.0), _Rshift(0.0) {};
    //! Construct potential from necessary constants
    ProteinMorse(double epsilon, double sigma, double Rshift): _epsilon(epsilon), _sigma(sigma), _Rshift(Rshift) {};

    //! Destructor
    virtual ~ProteinMorse() {};

    double computeEnergy(ProteinNode * A,  ProteinNode *B);

    void setScaling(double epsilon) { _epsilon = epsilon; };
    void setEpsilon(double epsilon) { _epsilon = epsilon; };
    void setSigma(double sigma) { _sigma = sigma; };
    void setShift(double Rshift) { _Rshift = Rshift; };  
    double getEpsilon() { return _epsilon; };
    double getSigma() { return _sigma; };
    double getShift() { return _Rshift; };
    
  private:
    // Members:
    double _epsilon;
    double _sigma;
    double _Rshift;
  }; 
  
}
#endif // _PROTEINMORSE_H_
