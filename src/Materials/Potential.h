// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file Potential.h

  \brief Base Class for a potential object.  Typical interaction potentials should be
  derived from this.

*/

#if !defined(__Potential_h__)
#define __Potential_h__

#include "voom.h"
#include "VoomMath.h"
#include "Node.h"

namespace voom
{

class Potential
{
 public:
  
  virtual void updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool f0, bool f1, bool f2) = 0;
  virtual double energy() const {return _W; };

  void ConsistencyTest(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, double eps = 1.0e-8, double tol = 1.0e-7);

  virtual void setScaling(double ) = 0;

  virtual double computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB) = 0;

 protected:
  double _W; // Energy

};

} //namespace voom

#endif //  !defined(__Potential_h__)
