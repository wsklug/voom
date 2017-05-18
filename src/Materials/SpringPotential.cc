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
#include "SpringPotential.h"

namespace voom {
        
  void SpringPotential::updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool fl0, bool fl1, bool fl2)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    if (fl0) {
      	// Spring energy
	_W = 0.5*_springK*pow(r - _restL, 2.0);
    }
    
    if (fl1) {
	// Spring force
	double factor = _springK*(r - _restL)/r;

	Vector3D ForceIncrement(0.0);
	ForceIncrement = factor*(nodeA->point() - nodeB->point());
	nodeA->updateForce(ForceIncrement);
	ForceIncrement = -ForceIncrement;
	nodeB->updateForce(ForceIncrement);
    }

    if (fl2) {
      // Not implemented 
      cout << "Stiffness calculations in spring potential not implemented " << endl;
      // exit(1);
    }

  } // SpringPotential::updateState

  double SpringPotential::computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB) {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    return _springK*(r - _restL);
  }

} // namespace voom
