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
// Revision 1.1  2005/10/22 19:11:22  klug
// Separated Element and StandardElement
//
//
//----------------------------------------------------------------------

#if !defined(__StandardElement_h__)
#define __StandardElement_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Element.h"

namespace voom
{

  //! StandardElement adds nodes and quadrature points to Element.
  /*! This class provides a generic template base for a standard
      finite element for different types of node, quadrature rule,
      material, and shape function.  Shape function and material
      objects are stored inside quadrature point objects.  This allows
      for convenient iteration through a container of quad points as
      is typically done for computation of element integrals.  The
      Quadrature_t object is to be used in initializing the weights
      and positions of the quad points.
   */
  template< class Node_t,
	    class Quadrature_t,
	    class Material_t,
	    class Shape_t >
  class StandardElement : public Element
  {
    
  public:
    
    //! default constructor
    StandardElement() {;}

    //! virtual destructor
    virtual ~StandardElement() {;}
    
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2) = 0;

    typedef typename std::vector<Node_t*> NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;

    //! Access the container of nodes
    virtual const NodeContainer& nodes() const { return _nodes; }

    struct QuadPointStruct 
    {
      double	 weight;
      Material_t material;
      Shape_t 	 shape;

      QuadPointStruct(double w, Material_t m, Shape_t s) 
	: weight(w), material(m), shape(s) {}
    };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
   		    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }

    QuadPointContainer & quadraturePoints() { return _quadPoints; }
    
  protected:
    NodeContainer _nodes;  
    QuadPointContainer _quadPoints;
  };
	
} // namespace voom

#endif // __StandardElement_h__
