// -*- C++ -*-
//----------------------------------------------------------------------
//
//                           Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file BrickElement.h

  \brief 3D brick FE

*/

#if !defined(__BrickElement_h__)
#define __BrickElement_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "MooneyRivlin.h"
#include "Quadrature.h"
#include "Shape.h"
#include "VoomMath.h"

namespace voom
{

  class BrickElement : public Element
  {

  public:

    typedef DeformationNode<3> Node; // nickname for mechanics nodes
   

    typedef std::vector< Node* > NodeContainer;
    typedef std::vector< Node* >::iterator NodeIterator;
    typedef std::vector< Node* >::const_iterator ConstNodeIterator;

    typedef MooneyRivlin MaterialType;

    //! virtual destructor
    virtual ~BrickElement() {;}

    //! Constructor

    BrickElement( const NodeContainer & Nodes,
		  MaterialType * material,
		  Quadrature<3> * quad,
		  Shape<3> * shape );

    struct QuadPointStruct 
    {
      double	 weight;
      MaterialType * material;
      Shape<3>::FunctionContainer shapeFunctions;
      Shape<3>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double w, MaterialType * m, Shape<3> * s, const NodeContainer & nds);
    };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef QuadPointContainer::iterator QuadPointIterator;
    typedef QuadPointContainer::const_iterator ConstQuadPointIterator;

    const QuadPointContainer & quadPoints() const {return _quadPoints;}

    const NodeContainer & Nodes() const {return _nodes;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }
		
    //! access volume
    double volume() const { return _volume; }

    //
    //   data
    //
  private:

    NodeContainer _nodes;

    Quadrature<3> * _quadrature;

    Shape<3> * _shape;

    QuadPointContainer _quadPoints;

    //! strain energy
    double _strainEnergy;

    //! volume corresponding to this element
    double _volume;

    blitz::Array< Vector3D, 1> _internalForce;

  };  // end of class
} // namespace voom

#endif // __BrickElement_h__
