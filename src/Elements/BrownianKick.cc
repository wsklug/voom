// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Created 2015/11/25 03:06:00  amit
// 
//
//
//----------------------------------------------------------------------


#include "VoomMath.h"
#include "BrownianKick.h"

namespace voom
{

  // Constructor
  BrownianKick::BrownianKick(const NodeContainer &defNodes, 
			     double Cd, double D, double dt )
    : _nodes(defNodes), _Cd(Cd), _D(D), _dt(dt) {
    // seed random number generator
    _rng.seed((unsigned int)time(0));
    // set the number of nodes
    _nodeCount = _nodes.size();
    _delta_xB.resize(_nodeCount);
  }

  void BrownianKick::updateKick(){
    for(int i=0; i < _nodeCount; i++){
      Vector3D currXi(_rng.random(),
				     _rng.random(),_rng.random());
      _delta_xB[i] = currXi*sqrt(_D*_dt);
    }
  }

  // Do mechanics on element; compute energy, forces, and/or stiffness.
  void BrownianKick::compute(bool f0, bool f1, bool f2) {
    
    if( f0 ) {
      _energy = 0;
      for(int i=0; i < _nodeCount; i++){
	double tempSum = 0;
	for(int j=0; j < 3; j++){
	  tempSum += _delta_xB[i][j]*_nodes[i]->getPoint(j);
	}
	_energy += -_Cd*tempSum;
      }
    }
    
    if( f1 ) {
      Vector3D f;
      for(int i=0; i < _nodeCount; i++){
      f = -_Cd * _delta_xB[i];
      _nodes[i]->updateForce(f);
      }
    }

    return;

  }
    
} // namespace voom
