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
  \file LennardJones.h
  \brief Interface for LennardJones potential in 3D cartesian space
*/

#ifndef _LENNARDJONES_H_
#define _LENNARDJONES_H_

#include "Potential.h"
#include "VoomMath.h"
#include "Node.h"

#include<set>

using namespace std;

namespace voom {
  
  class LennardJones : public Potential
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    LennardJones(): _epsilon(0.0), _sigma(0.0) {};
    //! Construct potential from necessary constants
    LennardJones(double epsilon, double sigma): _epsilon(epsilon), _sigma(sigma) {};

    //! Destructor
    virtual ~LennardJones() {}
    
    // Operators
    // General methods:
    //! Based on nodal postitions, calculates state of material (strain energy density and nodal forces)
    void updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool f0, bool f1, bool f2);

    double computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB);
  
    void setScaling(double epsilon) { _epsilon = epsilon; }; 
    void setEpsilon(double epsilon) { _epsilon = epsilon; };
    void setSigma(double sigma) { _sigma = sigma; };
    double getEpsilon() { return _epsilon; };
    double getSigma() { return _sigma; };
    
  private:
    
    // Members:
    double _epsilon;
    double _sigma;

  }; // LennardJones class
  
}
#endif // _LENNARDJONES_H_
