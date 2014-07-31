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
        
  void Potential::ConsistencyTest(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, double eps, double tol) {
    cout << endl << "Checking forces consistency in a potential class " << endl;
    // Clean nodal forces before consistency check	
    for (uint i = 0; i<3; i++) {
      	nodeA->setForce(i, 0.0);
	nodeB->setForce(i, 0.0);
    }
    // Compute current nodal forces
    updateState(nodeA, nodeB, false, true, false);
    
    Vector3D forceNumA(0.0), forceNumB(0.0);
    // Perturb center node position
    double W = 0.0;
    for(int i = 0; i < 3; i ++){
      nodeA->addPoint(i, eps);
      updateState(nodeA, nodeB, true, false, false);
      W = _W;
      nodeA->addPoint(i, -2.0*eps);
      updateState(nodeA, nodeB, true, false, false);
      W -= _W;
      nodeA->addPoint(i, eps);       // restore value
      forceNumA(i) = W/(2.0*eps);

      nodeB->addPoint(i, eps);
      updateState(nodeA, nodeB, true, false, false);
      W = _W;
      nodeB->addPoint(i, -2.0*eps);
      updateState(nodeA, nodeB, true, false, false);
      W -= _W;
      nodeB->addPoint(i, eps);       // restore value
      forceNumB(i) = W/(2.0*eps);
    }

    // Print results to screen
    cout << endl << "Analytical value of force vectors:" << endl;
    cout << nodeA->force() << endl;
    cout << nodeB->force() << endl;
    cout << "Numerical value of force vector:" << endl;
    cout << forceNumA << endl;
    cout << forceNumB << endl;
    double error =  sqrt(pow(tvmet::norm2(forceNumA - nodeA->force()), 2.0) + pow(tvmet::norm2(forceNumB - nodeB->force()), 2.0));
    cout << endl << "Numerical error:" << error << endl;
    double AbsTolerance = tol*sqrt(pow(tvmet::norm2(nodeA->force()), 2.0) + pow(tvmet::norm2(nodeB->force()), 2.0) );
    if (error < AbsTolerance ) {
      cout << "Potential consistency test PASSED. Error = " << error << " Tolerance = " << AbsTolerance << endl;
    }
    else {
      cout << "Potential consistency test FAILED. Error = " << error << " Tolerance = " << AbsTolerance << endl;
    }
    cout << endl;
    
  } //  Potential::ConsistencyTest


} //namespace voom

