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
  \file Morse.h
  \brief Interface for Morse potential in 3D cartesian space
*/

#ifndef _MORSE_H_
#define _MORSE_H_

#include "Potential.h"
#include "VoomMath.h"
#include "Node.h"

#include<set>

using namespace std;

namespace voom {
  
  class Morse : public Potential
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    Morse(): _epsilon(0.0), _sigma(0.0), _Rshift(0.0) {};
    //! Construct potential from necessary constants
    Morse(double epsilon, double sigma, double Rshift): _epsilon(epsilon), _sigma(sigma), _Rshift(Rshift) {};

    //! Destructor
    virtual ~Morse() {}
    
    // Operators
    // General methods:
    //! Based on nodal postitions, calculates state of material (strain energy density and nodal forces)
    void updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool f0, bool f1, bool f2);

    double computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB);

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

  }; // Morse class
  
}
#endif // _MORSE_H_
