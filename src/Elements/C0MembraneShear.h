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
  \file C0MembraneShear.h

  \brief C0MembraneShear is a concrete class derived from element, implementing
  a C0 finite element for finite deformation membrane structures with shear dof.

*/

#if !defined(__C0MembraneShear_h__)
#define __C0MembraneShear_h__

#include "voom.h"
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "Quadrature.h"
#include "Shape.h"
#include "EvansElastic_SkewedMin.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

using namespace std;

namespace voom
{

  /*!  Concrete class for a C0 finite deformation membrane finite element with shear dof */
 
  class C0MembraneShear : public Element 
  {
    public:

    typedef DeformationNode<3> Node;
    typedef vector< Node* > NodeContainer;
    typedef NodeContainer::iterator NodeIterator;

    //! virtual destructor
    virtual ~C0MembraneShear() {}

    C0MembraneShear(NodeContainer & Nodes,
		    ScalarFieldNode<3> * ShearNode,
		    ScalarFieldNode<3> * DirectionNode,
		    const double & AngleOffset,
		    const Quadrature<2> & Quad,
		    EvansElastic_SkewedMin * Mat,
		    Shape<2> * shape,
		    bool InsertShear = false,
		    bool InsertDirection = false);

    struct QuadPointStruct 
    {
      double	                 weight;
      EvansElastic_SkewedMin  *material;
      Shape<2>::FunctionContainer   shapeFunctions;
      Shape<2>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double wgt, EvansElastic_SkewedMin *mat, const Shape<2>::FunctionContainer & fun, const Shape<2>::DerivativeContainer & der):
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

    //! Return strain energy
    double strainEnergy() const { return _energy;}

    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    //! set the reference configuration in xy plane explicitly
    void SetRefConfiguration(double edgelen);



    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }

    const Tensor3D cauchyStress();
    const std::vector<double > matInvariants();

  private:
    // deformation nodes
    NodeContainer _nodes;

    // shear nodes
    ScalarFieldNode<3> * _shearNode;

    // direction nodes
    ScalarFieldNode<3> * _directionNode;

    // starting direction angle
    double _angleOffset;

    // quad points
    QuadPointContainer _quadPoints;

    // additional
    double _area, _volume;

  };  // end of class
} // namespace voom

// #include "C0MembraneShear.icc"

#endif // __C0MembraneShear_h__
