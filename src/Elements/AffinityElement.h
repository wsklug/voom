// -*- C++ -*-
//----------------------------------------------------------------------
//
//                              Mo Bai
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__AffinityElement_h__)
#define __AffinityElement_h__

#include <blitz/array.h>
#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "VoomMath.h"
#include "ShapeTri3.h"

namespace voom
{
  class AffinityElement 
    : public Element
  {
    
  public:
    
    typedef BrownianNode<2> Node_t;
    typedef ShapeTri3 Shape_t;

    typedef std::vector<Node_t*> NodeContainer;
    typedef NodeContainer::iterator NodeIterator;
    typedef NodeContainer::const_iterator ConstNodeIterator;

    //! virtual destructor
    virtual ~AffinityElement() {;}

  public:
    AffinityElement(const Shape_t & s, const NodeContainer & nodes ); 
    
    // Calculates strain tensor
    const Tensor2D & Strain() const;
    const Tensor2D & DisplacementGradient() const;
    double Rotation() const;
    double Area() const {return _area;}
    const Vector2D & Centroid();

    // Access the container of nodes
    const NodeContainer & nodes() const {return _nodes;}

    void compute();

    virtual void compute(bool f0, bool f1, bool f2) {;}


  private:

    NodeContainer _nodes;
    Shape_t _s;
    Vector2D _centroid;
    Tensor2D _strain;
    Tensor2D _dispGrad;
    double _rotation;
    double _area;
  };  // end of class
} // namespace voom

#endif // __AffinityElement_h__
