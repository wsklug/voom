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
    // return _epsilon*pow(1.0 - exp(-_sigma*(r - _Rshift)),2.0);
    return 0.5*_epsilon*( pow(1.0 - exp(-_sigma*(r - _Rshift)),2.0) - 1.0 );
   
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

  double ProteinMorse::computeForce(ProteinNode * A, ProteinNode *B){
    //TODO: Amit: Somebody needs to implement this if they want to use it
  }
  double ProteinMorse::computedWdEqPar(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // dW/d_Rshift
    // return 2.0*_epsilon*exp(_sigma*(_Rshift-r))*(-1.0 + exp(_sigma*(_Rshift-r)))*_sigma;
    return _epsilon*exp(_sigma*(_Rshift-r))*(-1.0 + exp(_sigma*(_Rshift-r)))*_sigma; // Should it be minus this expression ??
  }

  double ProteinMorse::computeddWddEqPar(ProteinNode * A,  ProteinNode *B)
  {
    double r = A->getDistance(B);
    // ddW/dd_Rshift
    // return 2.0*_epsilon*exp(_sigma*(_Rshift-r))*(-1.0 + 2.0*exp(_sigma*(_Rshift-r)))*pow(_sigma, 2.0);
    return _epsilon*exp(_sigma*(_Rshift-r))*(-1.0 + 2.0*exp(_sigma*(_Rshift-r)))*pow(_sigma, 2.0); // Should it be minus this expression ?
  }

} // namespace voom
