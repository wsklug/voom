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
#include "Morse.h"

namespace voom {
        
  void Morse::updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool fl0, bool fl1, bool fl2)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    if (fl0) {
      // Morse energy
      _W = _epsilon*(exp(-2*_sigma*(r-_Rshift))-2*exp(-_sigma*(r-_Rshift)));
    }
    
    if (fl1) {
      // Morse forces
      double factor = ((2.0*_sigma*exp(-(r-_Rshift)*_sigma)-2.0*_sigma*exp(-2.0*(r-_Rshift)*_sigma))
		       *_epsilon)/r;	
      Vector3D ForceIncrement(0.0);
      ForceIncrement = factor*(nodeA->point() - nodeB->point());
      nodeA->updateForce(ForceIncrement);
      ForceIncrement = -ForceIncrement;
      nodeB->updateForce(ForceIncrement);
    }

    if (fl2) {
      // Not implemented 
      cout << "Stiffness calculations in Morse potential not implemented " << endl;
      // exit(1);
    }

  } // Morse::updateState



  double Morse::computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    //return 2.0*_epsilon*_sigma*(1.0 - exp(-_sigma*(r - _Rshift)))*exp(-_sigma*(r - _Rshift));
    return (2.0*_sigma*exp(-(r-_Rshift)*_sigma)-2.0*_sigma*exp(-2.0*(r-_Rshift)*_sigma))
	    *_epsilon;
  } // Morse::computeTension


} // namespace voom
