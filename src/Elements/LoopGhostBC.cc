// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                    (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#include "LoopGhostBC.h"

namespace voom
{

  LoopGhostBC::LoopGhostBC(Node_t * N0, Node_t * N1, Node_t * N2, Node_t *N3) 
    : _N0(N0), _N1(N1), _N2(N2), _N3(N3) 
  {
    assert( _N0 != 0 );
    assert( _N1 != 0 );
    assert( _N2 != 0 );
    assert( _N3 != 0 );
  }
  
  void LoopGhostBC::predict() {
    const Vector3D & X0 = _N0->position();
    const Vector3D & X1 = _N1->position();
    const Vector3D & X2 = _N2->position();
    Vector3D X3;
    X3 = X0+X2-X1;
    _N3->setPosition(X3);

    const Vector3D & x0 = _N0->point();
    const Vector3D & x1 = _N1->point();
    const Vector3D & x2 = _N2->point();
    Vector3D x3;
    x3 = x0+x2-x1;
    _N3->setPoint(x3);
  }

  void LoopGhostBC::correct() {
    Vector3D f3 = _N3->force();
    _N0->updateForce(f3);
    _N2->updateForce(f3);
    f3 = -f3;
    _N1->updateForce(f3);

    for(int i=0; i<_N3->dof(); i++) _N3->setForce(i,0.0);
  }

} // namespace voom
