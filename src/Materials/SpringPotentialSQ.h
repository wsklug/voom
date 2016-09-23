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
  \file SpringPotentialSQ.h
  \brief Interface for spring potential in 3D cartesian space
*/

#ifndef _SPRINGPOTENTIALSQ_H_
#define _SPRINGPOTENTIALSQ_H_

#include "Potential.h"
#include "VoomMath.h"
#include "Node.h"

#include<set>

using namespace std;

namespace voom {
  
  class SpringPotentialSQ : public Potential
  {
    
  public:
    
    // Constructors/destructors:
    //! Default constructor
    SpringPotentialSQ(): _springK(0.0), _restL(0.0) {};
    //! Construct potential from necessary constants
    SpringPotentialSQ(double springK, double restL): _springK(springK), _restL(restL) {};

    //! Destructor
    virtual ~SpringPotentialSQ() {}
    
    // Operators
    // General methods:
    //! Based on nodal postitions, calculates state of material (strain energy and nodal forces)
    void updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool f0, bool f1, bool f2);

    double computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB);

    void setScaling(double springK) { _springK = springK; };   
    void setSpringK(double springK) { _springK = springK; };
    void setRestL(double restL) { _restL = restL; };
    double getSpringK() { return _springK; };
    double getRestL() { return _restL; };
    
  private:
    
    // Members:
    double _springK;
    double _restL;

  }; // SpringPotentialSQ class
  
}
#endif // _SPRINGPOTENTIALSQ_H_
