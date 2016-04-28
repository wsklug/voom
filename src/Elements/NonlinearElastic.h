// -*- C++ -*-
//----------------------------------------------------------------------
//
//                       Melissa M. Gibbons
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file NonlinearElastic.h

  \brief Nonlinear Elastic is a concrete class derived from Element, implementing
  a nonlinear elastic finite element for finite deformation of 3D structures.

*/

#if !defined(__NonlinearElastic_h__)
#define __NonlinearElastic_h__

#include <blitz/array.h>
 
#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Concrete class for a finite deformation 3D finite element.
  */
  template<class Quadrature_t, class Material_t, class Shape_t>
  class NonlinearElastic 
    : public Element
  {
    
  public:
    
    typedef DeformationNode<3> Node_t;

    typedef typename std::vector<Node_t*> NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;

    //! virtual destructor
    virtual ~NonlinearElastic() {;}

    //! Constructor; from quadrature and material objects and node container
    NonlinearElastic( const Quadrature_t & quad,
		      const Material_t & mat,
		      const NodeContainer & nodes
		      );

  public:
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    // helper functions for post-processing

    //! Calculates von Mises stress in an element
    double CalcMisesStress();
    //! Calculates the Cauchy stress tensor in an element
    Tensor3D CalcCauchyStress();
    //! Calculates the principal strain vector in an element
    Vector3D CalcPrincipalStrains();
  
    void setYoungsModulus(double value);
    // Accessor functions

    //! Returns strain energy stored in an element
    double strainEnergy() const { return _strainEnergy; }
    //! Returns element volume
    double elVolume() const { return _volume; }

    //! Access the container of nodes
    const NodeContainer & nodes() { return _nodes; }
    
    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();


    struct QuadPointStruct 
    {
      double	 weight;
      typename Shape_t::FunctionContainer shapeFunctions;
      typename Shape_t::DerivativeContainer shapeDerivatives;
      Material_t material;
      QuadPointStruct(double w, const Material_t & m, const Shape_t & s, const NodeContainer & nds);
    };

    typedef typename std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }
    QuadPointContainer & quadraturePoints() { return _quadPoints; }

    //
    //   data
    //
  private:
    //! Forces on the element
    blitz::Array< Vector3D, 1> _internalForce;

    //! Strain energy
    double _strainEnergy;

    //! Element volume
    double _volume;

    //! Node container
    NodeContainer _nodes;

    //! Container of quadrature points
    QuadPointContainer _quadPoints;

  };  // end of class
} // namespace voom

#include "NonlinearElastic.icc"

#endif // __NonlinearElastic_h__
