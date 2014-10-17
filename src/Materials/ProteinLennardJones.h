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

#ifndef _PROTEINLENNARDJONES_H_
#define _PROTEINLENNARDJONES_H_

#include "ProteinPotential.h"

using namespace std;

namespace voom {
  
  class ProteinLennardJones: public ProteinPotential
  {
  public:
    // Constructors/destructors:
    //! Default constructor
    ProteinLennardJones(): _epsilon(0.0), _sigma(0.0) {};
    //! Construct potential from necessary constants
    ProteinLennardJones(double epsilon, double sigma): _epsilon(epsilon), _sigma(sigma) {};

    //! Destructor
    virtual ~ProteinLennardJones() {};

    double computeEnergy(ProteinNode * A,  ProteinNode *B);
    double computedWdEqPar(ProteinNode * A,  ProteinNode *B);
    double computeddWddEqPar(ProteinNode * A,  ProteinNode *B);

    void setEquilibriumParam(double a){ _sigma = a; };
    double getEquilibriumParam(){ return _sigma; };
    void setEquilibriumR(double R){ _sigma = R/pow(2.0, 1.0/6.0); };
    double getEquilibriumR(){ return (_sigma*pow(2.0, 1.0/6.0)); };

    void setScaling(double epsilon) { _epsilon = epsilon; };
    void setEpsilon(double epsilon) { _epsilon = epsilon; };
    void setSigma(double sigma) { _sigma = sigma; }; 
    double getEpsilon() { return _epsilon; };
    double getSigma() { return _sigma; };
    
  private:
    // Members:
    double _epsilon;
    double _sigma;
  }; 
  
}
#endif // _PROTEINLENNARDJONES_H_
