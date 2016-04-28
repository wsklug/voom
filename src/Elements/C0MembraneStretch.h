// -*- C++ -*-
//----------------------------------------------------------------------
//
//                     William S. Klug, Feng Feng
//                     Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file C0MembraneStretch.h

  \brief C0MembraneStretch is a concrete class derived from element, implementing
  a C0 finite element for finite deformation membrane structures with stretch dof.

*/

#if !defined(__C0MembraneStretch_h__)
#define __C0MembraneStretch_h__

#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "Quadrature.h"
#include "Shape.h"
#include "EvansElastic_Stretch.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

using namespace std;

namespace voom
{

  /*!  Concrete class for a C0 finite deformation membrane finite element with stretch dof */
 
  class C0MembraneStretch : public Element 
  {
    public:

    typedef DeformationNode<3> Node;
    typedef vector< Node* > NodeContainer;
    typedef NodeContainer::iterator NodeIterator;

    //! virtual destructor
    virtual ~C0MembraneStretch() {}

    C0MembraneStretch(NodeContainer & Nodes,
		      ScalarFieldNode<3> * StretchNode,
		      ScalarFieldNode<3> * DirectionNode,
		      ScalarFieldNode<3> * phiNode,
		      const double & AngleOffset,
		      const Quadrature<2> & Quad,
		      EvansElastic_Stretch * Mat,
		      Shape<2> * shape,
		      bool InsertStretch = false,
		      bool InsertDirection = false,
		      bool InsertPhi = false);

    struct QuadPointStruct 
    {
      double	               weight;
      EvansElastic_Stretch  *material;
      Shape<2>::FunctionContainer   shapeFunctions;
      Shape<2>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double wgt, EvansElastic_Stretch *mat, const Shape<2>::FunctionContainer & fun, const Shape<2>::DerivativeContainer & der):
      weight(wgt), material(mat)      
      {
	shapeFunctions = fun;
	shapeDerivatives = der;
      }
    };

    typedef vector<QuadPointStruct> QuadPointContainer;
    typedef QuadPointContainer::iterator QuadPointIterator;
    typedef QuadPointContainer::const_iterator ConstQuadPointIterator;


    const NodeContainer & nodes() const {return _nodes;}
    const QuadPointContainer & quadPoints() const {return _quadPoints;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    //! Calculate and return the two surface deformation invariants  
    void invariants(double& I1, double& J);

    //! Push forward operator
    Vector3D PushForwardOperator(Vector3D & Nbar);

    //! Return strain energy
    double strainEnergy() const { return _energy;}
    double conformationalEnergy() const { return _confEnergy;}

    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    //! set the reference configuration in xy plane explicitly
    void SetRefConfiguration(double edgelen);

    //! reset the angle offset
    void setAngleOffset(double new_AngleOffset) { _angleOffset = new_AngleOffset;}

    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }

    const Tensor3D cauchyStress();
    const vector<double > matInvariants();

  private:
    // deformation nodes
    NodeContainer _nodes;

    // stretch node
    ScalarFieldNode<3> * _stretchNode;

    // direction node
    ScalarFieldNode<3> * _directionNode;

    // soft mode direction node
    ScalarFieldNode<3> * _phiNode;

    // starting direction angle
    double _angleOffset;

    // quad points
    QuadPointContainer _quadPoints;

    // additional
    double _area, _volume, _confEnergy;

  };  // end of class
} // namespace voom

// #include "C0MembraneStretch.icc"

#endif // __C0MembraneStretch_h__
