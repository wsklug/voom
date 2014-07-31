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
    return 4.0*_epsilon*(pow(_sigma/r, 12.0) - pow(_sigma/r, 6.0)); 
  }

} // namespace voom
