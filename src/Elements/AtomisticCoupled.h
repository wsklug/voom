// -*- C++ -*-
//----------------------------------------------------------------------
//
//                       Melissa M. Gibbons
//                University of California Los Angeles
//                   (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file AtomisticCoupled.h

  \brief Atomistic Coupled is a concrete class derived from element, implementing
  an atomistically coupled finite element for finite deformation of 3D structures,
  wherein the forces are calculated at the atomic level, using MD codes.

*/

#if !defined(__AtomisticCoupled_h__)
#define __AtomisticCoupled_h__

#include <blitz/array.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"
#include "Math.h"

namespace voom
{

  /*!  Concrete class for a finite deformation 3D finite element, where
    forces are calculated at the atomic level.
  */
  template<class Shape_t>
  class AtomisticCoupled 
    : public Element
  {
    
  public:
    // typedefs
    typedef tvmet::Vector<double, 3> Vector3D;
    typedef tvmet::Matrix<double, 3, 3> Matrix3D;

    typedef DeformationNode<3> Node_t;

    typedef typename std::vector<Node_t*> NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;

    //typedef typename std::vector< AtomCoord > AtomInContainer;
    typedef typename std::vector< DeformationNode<3>* > AtomContainer;
    typedef typename AtomContainer::iterator AtomIterator;
    typedef typename AtomContainer::const_iterator ConstAtomIterator;

    typedef typename tvmet::Vector<double,3> Point;

    //! virtual destructor
    virtual ~AtomisticCoupled() {;}

    //! Constructor
    AtomisticCoupled(const NodeContainer & nodes, const AtomContainer & atoms);

  public:
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);

    // Based on new nodal positions, update and set new atom positions
    void recalcAtomGlobPos();
    
    // helper functions for post-processing
    //Vector3D CalcPrincipalStrains();

    // Accessor functions
    double strainEnergy() const { return _strainEnergy; }

    //! Access the container of nodes
    const NodeContainer & nodes() { return _nodes; }
    
    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();

    //! Access the container of atoms
    const AtomContainer & atoms() const { return _atoms; }
    AtomContainer & atoms() { return _atoms; }

    
    //
    //   data
    //
  private:
    //! forces on the element
    blitz::Array< Vector3D, 1> _internalForce;

    //! strain energy
    double _strainEnergy;

    //! nodes contained in the element
    NodeContainer _nodes;

    // Atoms contained in the element
    AtomContainer _atoms;

    // element shape functions
    blitz::Array<double, 2> _shapeFunctions;

  };  // end of class
} // namespace voom

#include "AtomisticCoupled.icc"

#endif // __AtomisticCoupled_h__
