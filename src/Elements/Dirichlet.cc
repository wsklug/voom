// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#include "Dirichlet.h"

namespace voom
{

  void Dirichlet::predict() {
    _node->setPoint(_dof,_value);
    return;
  }

  void Dirichlet::correct() {
    _node->setForce(_dof, 0.0);
    return;
  }

  void Dirichlet::updateValue() {
    _value = _node->getPoint(_dof);
  }

} // namespace voom
