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
#include "LennardJones.h"

namespace voom {
        
  void LennardJones::updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool fl0, bool fl1, bool fl2)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    if (fl0) {
      	// Lennard Jones energy
      	_W = 4.0*_epsilon*(pow(_sigma/r, 12.0) - pow(_sigma/r, 6.0)); 
    }
    
    if (fl1) {
	// Lennard Jones forces
	double factor = (4.0*_epsilon/pow(r, 2.0))*(6.0*pow(_sigma/r, 6.0) - 12.0*pow(_sigma/r, 12.0));
	
        Vector3D ForceIncrement(0.0);
	ForceIncrement = factor*(nodeA->point() - nodeB->point());
	nodeA->updateForce(ForceIncrement);
	ForceIncrement = -ForceIncrement;
	nodeB->updateForce(ForceIncrement);
    }

    if (fl2) {
      // Not implemented 
      cout << "Stiffness calculations in LennardJones potential not implemented " << endl;
      // exit(1);
    }

  } // LennardJones::updateState


  double LennardJones::computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    return (4.0*_epsilon/r)*(6.0*pow(_sigma/r, 6.0) - 12.0*pow(_sigma/r, 12.0));
  } // LennardJones::computeTension



} // namespace voom
