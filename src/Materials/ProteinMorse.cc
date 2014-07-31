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
#include "ProteinMorse.h"

namespace voom {
        
  double ProteinMorse::computeEnergy(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // Morse potential type energy
    return _epsilon*pow(1.0 - exp(-_sigma*(r - _Rshift)),2.0);
   
    // if (fl1) {
    //   // Morse forces
    //   double factor = 2.0*_epsilon*_sigma*(1.0 - exp(-_sigma*(r - _Rshift)))*exp(-_sigma*(r - _Rshift))/r;
	
    //     Vector3D ForceIncrement(0.0);
    // 	ForceIncrement = factor*(nodeA->point() - nodeB->point());
    // 	nodeA->updateForce(ForceIncrement);
    // 	ForceIncrement = -ForceIncrement;
    // 	nodeB->updateForce(ForceIncrement);
    // }

    // if (fl2) {
    //   // Not implemented 
    //   cout << "Stiffness calculations in LennardJonesFT potential not implemented " << endl;
    //   // exit(1);
    // }

  } // Morse::updateState



} // namespace voom
