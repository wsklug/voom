// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
#include "S2.h"

namespace voom
{
  void S2::predict() 
  {
    Vector3D x;
    x = _node->point();
    double R =  norm2(x-_xc);
    if( R==0 ) {
      std::cout << "S2::predict(): R==0.  Exiting." << std::endl;
      exit(0);
    }
    
    if( R != _R ) {
      // Node not on the sphere of radius R centered at xc.  Project it.
      double dR = R - _R;
      Vector3D dx;
      dx = dR*(x-_xc)/R;
      x -= dx;
      _node->setPoint(x);
    } 
    return;
    }
  
  void S2::correct() 
  {
    // add reaction to balance force in the radial direction
    Vector3D n;
    n = _node->point() - _xc;
    n /=  norm2(n);
    _f = _node->force();
    double fn =  dot(_f,n);
    _f = -fn*n;
    _node->updateForce( _f );
    return;
  }

  double S2::getForce(int a, int i) const 
  {
    assert( a == 0 );
    assert( i >= 0 && i < 3 );
    return _f(i);
  }

}; // namespace voom

