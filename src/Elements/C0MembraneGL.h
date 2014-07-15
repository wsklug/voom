// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file C0MembraneGL.h

  \brief C0MembraneGL is a concrete class derived from element, implementing
  a C0 finite element for finite deformation membrane structures.

*/

#if !defined(__C0MembraneGL_h__)
#define __C0MembraneGL_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "GLElastic.h"
#include "Quadrature.h"
#include "Shape.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Element coupling membrane mechanics to a Ginzburg-Landau order-parameter scalar field.
  */
  class C0MembraneGL : public Element
  {

  public:

    typedef DeformationNode<3> DefNode; // nickname for mechanics nodes
    typedef ScalarFieldNode<3> GLNode;  // nickname for order-parameter nodes

    typedef std::vector< DefNode* > DefNodeContainer;
    typedef std::vector< DefNode* >::iterator DefNodeIterator;
    typedef std::vector< DefNode* >::const_iterator ConstDefNodeIterator;

    typedef std::vector< GLNode* > GLNodeContainer;
    typedef std::vector< GLNode* >::iterator GLNodeIterator;
    typedef std::vector< GLNode* >::const_iterator ConstGLNodeIterator;

    typedef GLElastic MaterialType;

    //! virtual destructor
    virtual ~C0MembraneGL() {;}

    //! Constructor

    C0MembraneGL( const DefNodeContainer & defNodes,
		  const GLNodeContainer  & glNodes ,
		  const MaterialType * material,
		  Quadrature<2> * quad,
		  Shape<2> * shape );

    struct QuadPointStruct 
    {
      double	 weight;
      MaterialType material;
      Shape<2>::FunctionContainer shapeFunctions;
      Shape<2>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double w, const MaterialType * m, Shape<2> * s);
    };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef QuadPointContainer::iterator QuadPointIterator;
    typedef QuadPointContainer::const_iterator ConstQuadPointIterator;

    const QuadPointContainer & quadPoints() const {return _quadPoints;}

    const DefNodeContainer & defNodes() const {return _defNodes;}

    const GLNodeContainer & glNodes() const {return _glNodes;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }

    //! Accessor for element work done by pressure
    double work() const { return _work; }
		
    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }

    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    //
    //   data
    //
  private:

    DefNodeContainer _defNodes;

    GLNodeContainer _glNodes;

    Quadrature<2> * _quadrature;

    Shape<2> * _shape;

    QuadPointContainer _quadPoints;

    //! strain energy
    double _strainEnergy;

    //! work done by pressure
    double _work;

    //! volume corresponding to this element
    double _volume;

    //! area enclosed by the current element
    double _area;

  };  // end of class
} // namespace voom

#endif // __C0MembraneGL_h__
