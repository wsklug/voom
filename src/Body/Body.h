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
// Revision 1.21  2005/10/22 19:15:40  klug
// Separated Node and NodeBase.  GCC 4 compatibility fixes.
//
// Revision 1.20  2005/10/21 01:35:11  klug
// Added pushBack method.
//
// Revision 1.19  2005/08/22 22:08:41  klug
// Assembly shifted from elements to nodes.  Bodies now assemble energy.
//
// Revision 1.18  2005/08/20 01:28:31  klug
// acinclude.m4
//
// Revision 1.17  2005/05/25 02:13:56  klug
// Removed Constraint.h Potential.h includes.
//
//

/*! 
  \file Body.h

  \brief Body is virtual base class specifying the interface for
  classes implementing the concept of a body on which to do mechanics.
  A body is composed of elements which are all of the same type; for
  each new element type, a new concrete class should be derived from
  Body.

*/

#if !defined(__Body_h__)
#define __Body_h__

#include<blitz/array.h>
#include<vector>
#include "NodeBase.h"
#include "Element.h"
#include "Constraint.h"

namespace voom
{

  /*!  Virtual base class for a body composed of elements all of the
    same type; classes for specific types of elements should be
    derived from Body
  */
  class Body
  {
    
  public:

    typedef std::vector<Element*> ElementContainer;
    typedef ElementContainer::iterator ElementIterator;
    typedef ElementContainer::const_iterator ConstElementIterator;

    typedef std::vector<Constraint*> ConstraintContainer;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef ConstraintContainer::const_iterator ConstConstraintIterator;

    typedef std::vector<NodeBase*> NodeContainer;
    typedef NodeContainer::iterator NodeIterator;
    typedef NodeContainer::const_iterator ConstNodeIterator;
    
    //! Default Constructor
    Body() {_output=paraview; _energy=0.0;}
    
    //! Default Destructor
    virtual ~Body() {};
    
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) = 0;
    
    double energy() const { return _energy; }

    //! Get the number of degrees of freedom in the body
    int dof() const {return _dof;} 
    
    ElementContainer & elements() {return _elements;}

    ConstraintContainer & constraints() {return _constraints;}

    const NodeContainer & nodes() const {return _nodes;}
    NodeContainer & nodes() {return _nodes;}

    //! Add an element to the list
    virtual void addElement( Element * e ) { _elements.push_back( e ); }

    //! Add a node to the list
    virtual void addNode( NodeBase * n ) { 
      _nodes.push_back( n ); 
      _dof += n->dof();
    }

    //! Add a constraint to the list
    virtual void addConstraint( Constraint * e ) { _constraints.push_back( e ); }

    //! Query an element's activity status
	virtual bool active(int e) {};

    //! Mark an element as active so it will be computed
    virtual void activate(int e) {};

    //! Mark an element as inactive so it will not be computed
    virtual void deactivate(int e) {};

    void setOutput( OutputSelection o ) { _output = o; }

    //! Print results
    virtual void print(std::string name) const {
      if(_output==paraview) printParaview(name);
    };

    //! Print input file for Paraview
    virtual void printParaview(std::string name) const = 0;

    virtual void checkConsistency(bool verbose=false);

  protected:

    int _id;

    int _dof;

    OutputSelection _output;

    ElementContainer _elements;

    NodeContainer _nodes;

    ConstraintContainer _constraints;

    double _energy;
    
  };

} // namespace voom
#endif // __Body_h__
