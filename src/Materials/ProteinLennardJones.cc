// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "ProteinLennardJones.h"

namespace voom {
        
  double ProteinLennardJones::computeEnergy(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // Lennard Jones energy
    // return 4.0*_epsilon*(pow(_sigma/r, 12.0) - pow(_sigma/r, 6.0)); 
    return 2.0*_epsilon*(pow(_sigma/r, 12.0) - pow(_sigma/r, 6.0)); // only return half energy - the rest comes from counting the other particle
  }

  double ProteinLennardJones::computeForce(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // Lennard Jones forces
    double factor = (4.0*_epsilon/pow(r, 2.0))*(6.0*pow(_sigma/r, 6.0) - 12.0*pow(_sigma/r, 12.0));
    
    Vector3D ForceIncrement(0.0);
    ForceIncrement = factor*(A->getHostPosition() - B->getHostPosition());
    A->getHost()->updateForce(ForceIncrement);
    ForceIncrement = -ForceIncrement;
    B->getHost()->updateForce(ForceIncrement);
  }

  double ProteinLennardJones::computedWdEqPar(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // Lennard Jones first derivative with respect to sigma
    // return 4.0*_epsilon*(12.0*pow(_sigma/r, 11.0) - 6.0*pow(_sigma/r, 5.0))/r; 
    return 2.0*_epsilon*(12.0*pow(_sigma/r, 11.0) - 6.0*pow(_sigma/r, 5.0))/r; 
  }

  double ProteinLennardJones::computeddWddEqPar(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // Lennard Jones energy
    // return 4.0*_epsilon*(132.0*pow(_sigma/r, 10.0) - 30.0*pow(_sigma/r, 4.0))/(r*r); 
    return 2.0*_epsilon*(132.0*pow(_sigma/r, 10.0) - 30.0*pow(_sigma/r, 4.0))/(r*r); 
  }

} // namespace voom
