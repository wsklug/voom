// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------


#if !defined(__MembraneGLImplicitMass_h__)
#define __MembraneGLImplicitMass_h__

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

  /*!  Mass Element coupling membrane mechanics to a Ginzburg-Landau order-parameter scalar field.
  */

  class MembraneGLImplicitMass : public Element
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

    //! virtual destructor
    virtual ~MembraneGLImplicitMass() {;}

    //! Constructor

    MembraneGLImplicitMass( const DefNodeContainer & defNodes,
			    const GLNodeContainer  & glNodes ,
			    Quadrature<2> * quad,
			    Shape<2> * shape,
			    double Mx0, double Mx1, double Meta0 );

    struct QuadPointStruct 
    {
      double	 weight;
      Shape<2>::FunctionContainer shapeFunctions;
      Shape<2>::DerivativeContainer shapeDerivatives;
      QuadPointStruct(double w, Shape<2> * s);
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

    //! step in time; save current nodal values as previous, and
    //! recompute masses

    void step(double dt);

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

    double _Mx0, _Mx1, _Meta0;

    //! mass matrices
    Array2D _Mx, _Meta;

    //! previous nodal positions
    std::vector< Vector3D > _x_prev;

    //! previous nodal densities
    std::vector< double > _eta_prev;

  };  // end of class
} // namespace voom

#endif // __MembraneGLImplicitMass_h__
