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
// Revision 1.13  2005/10/22 19:11:23  klug
// Separated Element and StandardElement
//
// Revision 1.12  2005/10/21 01:07:21  klug
// Cleaned and removed old stuff related to storing forces and stiffness
// in elements.  Note: stiffness part needs to be fixed in checkConsistency()
// and checkRank() is totally broken.
//
// Revision 1.11  2005/08/22 22:18:30  klug
// Assembly shifted from elements to nodes.  Element class renamed
// StandardElement, ElementBase class renamed Element, and verification
// tests moved to Element class and their definitions moved from Element.icc
// to Element.cc which is now compiled into libElement.a.
//
// Revision 1.10  2005/05/23 17:25:08  klug
// Added index map; renamed strainEnergy as energy.
//
// Revision 1.9  2005/04/11 16:24:12  klug
// New Element base classes with forces and stiffness stored locally.
//
//----------------------------------------------------------------------

/*! 
  \file Element.h

  \brief Class for a finite element.

*/

#if !defined(__Element_h__)
#define __Element_h__

#include <blitz/array.h>
#include <vector>
#include "voom.h"
#include "NodeBase.h"

namespace voom
{

  //! Virtual base class for a finite element.
  /*!  A generic element provides the capability for calculation of
    the strain energy for the element domain, and the first and second
    deriviatives of that energy with respect to the degrees of freedom
    of nodes attached to the element.  These calculations are
    performed by the method <tt>compute(bool f0, bool f1, bool
    f2)</tt> which has the behavior
    <ul>
    <li> If <tt> f0==true </tt> compute energy.
    <li> If <tt> f1==true </tt> compute first derivatives of energy (forces).
    <li> If <tt> f2==true </tt> compute second derivatives of energy (stiffness).
    </ul>

    Element also provides to testing methods,
    <tt>checkConsistency()</tt> which tests the numerical consistency
    of element derivatives, and <tt>checkRank()</tt> which tests the
    rank of the stiffness.
    
  */
  class Element
  {
  public:

    typedef std::vector<NodeBase*> BaseNodeContainer;
    typedef std::vector<NodeBase*>::iterator BaseNodeIterator;
    typedef std::vector<NodeBase*>::const_iterator ConstBaseNodeIterator;

    virtual ~Element() {}
    
    //! Accessor for element strain energy
    double energy() const { return _energy; } 

    //! Accessor for container of base nodes
    const BaseNodeContainer & baseNodes() const {return _baseNodes;}

    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2) = 0;
    
    virtual bool checkConsistency();
		
    //! Rank of the element stiffness matrix
    virtual bool checkRank(const int rank);

    virtual const Tensor3D cauchyStress()
    {
      Tensor3D A(0.0); 
      return A;
    }; 

    virtual const std::vector<double > matInvariants()
    {
      std::vector<double > I(2, 0.0); 
      return I;
    };

  protected:
    double _energy;

    BaseNodeContainer _baseNodes;

  };

	
} // namespace voom

#endif // __Element_h__
