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
#include "LennardJonesFT.h"

namespace voom {
        
  void LennardJonesFT::updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool fl0, bool fl1, bool fl2)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    if (fl0) {
      // Lennard JonesFT energy
      _W = _epsilon*(pow(_sigma/(r-_Rshift), 4.0) - pow(_sigma/(r-_Rshift), 2.0)); 
    }
    
    if (fl1) {
	// Lennard JonesFT forces
      double factor = (_epsilon/(r-_Rshift))*(2.0*pow(_sigma/(r-_Rshift), 2.0) - 4.0*pow(_sigma/(r-_Rshift), 4.0))/r;
	
        Vector3D ForceIncrement(0.0);
	ForceIncrement = factor*(nodeA->point() - nodeB->point());
	nodeA->updateForce(ForceIncrement);
	ForceIncrement = -ForceIncrement;
	nodeB->updateForce(ForceIncrement);
    }

    if (fl2) {
      // Not implemented 
      cout << "Stiffness calculations in LennardJonesFT potential not implemented " << endl;
      // exit(1);
    }

  } // LennardJonesFT::updateState


  double LennardJonesFT::computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB)
  {
    double r = tvmet::norm2( nodeA->point() - nodeB->point() );
    return (_epsilon/(r-_Rshift))*(2.0*pow(_sigma/(r-_Rshift), 2.0) - 4.0*pow(_sigma/(r-_Rshift), 4.0));
  } // LennardJones::computeTension


} // namespace voom
