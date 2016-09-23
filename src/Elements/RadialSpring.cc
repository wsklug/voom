// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2016 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Created 2016/08/04 19:47:00  amit
// 
//
//
//----------------------------------------------------------------------


#include "VoomMath.h"
#include "RadialSpring.h"

namespace voom
{

  // Constructor
  RadialSpring::RadialSpring(const NodeContainer &defNodes, 
			     double k, double R)
    : _nodes(defNodes), _k(k), _R(R) {

    // set the number of nodes
    _nodeCount = _nodes.size();
  }
  
  //Destructor
  RadialSpring::~RadialSpring(){}

  // Do mechanics on element; compute energy, forces, and/or stiffness.
  void RadialSpring::compute(bool f0, bool f1, bool f2) {

    //Calculate centre of sphere as average of position vectors of all nodes.
    Vector3D Xavg(0.0);
    for ( int i = 0; i < _nodeCount; i++){
      Xavg += _nodes[i]->point();
    }
    Xavg /= _nodeCount;

    double r;
    
    if( f0 ) {
      _energy = 0;
      for(int i=0; i < _nodeCount; i++){
	r = tvmet::norm2(_nodes[i]->point() - Xavg);
	_energy += (0.5*_k*(r - _R)*(r - _R));
      }
    }
    
    if( f1 ) {
      Vector3D f;
      for(int i = 0; i < _nodeCount; i++){
	r = tvmet::norm2(_nodes[i]->point() - Xavg);
	f = (_k*(r - _R)/r)*(_nodes[i]->point() - Xavg);
      _nodes[i]->updateForce(f);
      }
    }

    return;

  }
    
} // namespace voom
