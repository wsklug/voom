// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Dirichlet_h__)
#define __Dirichlet_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Constraint.h"

namespace voom
{

  //! Predictor-corrector Constraint class for a Derichlet boundary condition
  class Dirichlet : public Constraint
  {
  public:

    //! Constructor
    Dirichlet( NodeBase *  node, int i ) : _node(node), _dof(i) 
    { updateValue(); }

    //! Constructor
    Dirichlet( NodeBase *  node, int i, double value ) 
      : _node(node), _dof(i), _value(value) {}

    //! Adjust DOF to meet constraint; assume the initial 
    virtual void predict();

    //! Correct out of balance force conjugate to DOF
    virtual void correct();

    //! Access and store new value from the node
    void updateValue();

  private:

    //! The node to be constrained
    NodeBase * _node;

    //! The dof to be fixed
    int _dof;

    //! The value of the fixed dof
    double _value;
  };

} // namespace voom

#endif // __Dirichlet_h__
