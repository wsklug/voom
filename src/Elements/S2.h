// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__S2_h__)
#define __S2_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Constraint.h"

namespace voom
{

  class S2 : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<3> DefNode;

    S2(  DefNode * node, double R, const Vector3D & xc )
      : _node(node), _R(R), _xc(xc), _f(0.0)
    {}


    S2(  DefNode * node, double R=1.0 )
      : _node(node), _R(R), _xc(0.0), _f(0.0)
    {}

  

    virtual void predict(); 

    virtual void correct();

    double getForce(int i) const {return _f(i);}

    const Vector3D & force() const {return _f;}

    double getForce(int a, int i) const;

  private:

    double _R;
    Vector3D _xc;
    Vector3D _f;
    DefNode * _node;

  };



}; // namespace voom

#endif // __S2_h__
