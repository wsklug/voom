// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug & Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/10/22 19:13:37  klug
// Separated Node and NodeBase.  GCC 4 compatibility fixes.
//
//----------------------------------------------------------------------


#ifndef _NodeBase_
#define _NodeBase_

#ifdef _X
#undef _X
#endif
#ifdef _A
#undef _A
#endif
#ifdef _B
#undef _B
#endif

 
#include <vector>
#include <fstream>
#include "voom.h"

namespace voom
{
  class NodeBase 
  {
  public:
    typedef std::vector<int> DofIndexMap;
    typedef std::vector<NodeBase*> NeighborContainer;
    typedef std::vector<NodeBase*>::iterator NeighborIterator;
    typedef std::vector<NodeBase*>::const_iterator ConstNeighborIterator;

    NodeBase(int id, const DofIndexMap & index) : _id(id), _index(index) {;}

    virtual ~NodeBase() {}

    virtual int dof() const = 0;
    virtual double getPoint(int i) const = 0;
    virtual void setPoint(int i, double x) = 0;
    virtual void addPoint(int i, double dx) = 0;    

    virtual double getForce(int i) const = 0;
    virtual void setForce(int i, double f) = 0;
    virtual void addForce(int i, double f) = 0;    

    virtual double getStiffness(const NodeBase* b, int kb, int ia) const {
      return 0;
    }
    virtual void setStiffness(const NodeBase* b, int kb, int ia, double k) {;}
    virtual void addStiffness(const NodeBase* b, int kb, int ia, double k) {;}    
    virtual double getStiffness(int b, int kb, int ia) const {return 0;};
    virtual void setStiffness(int b, int kb, int ia, double k) {;};
    virtual void addStiffness(int b, int kb, int ia, double k) {;};    

    virtual double getStiffness(int ia) const { return 0; }
    virtual void setStiffness(int ia, double k) {};
    virtual void addStiffness(int ia, double k) {};    

    virtual double getMass() const {return 0;}

    virtual void setMass(double m) {};
    
    virtual void addMass(double m) {};
    
    
    //! assign id
    virtual void setId( const int i ) { _id = i; }

    //! access id
    virtual int id() const { return _id; }

    //! Mutator for global DOF index map
    void setIndex(const DofIndexMap & index) {_index = index;}

    //! Accessor for global DOF index map
    const DofIndexMap & index() const {return _index;}

    const NeighborContainer & neighbors() const {return _neighbors;}

  protected:
    int _id;
    DofIndexMap _index;
    NeighborContainer _neighbors;
  };

}; // namespace voom

#endif // _NodeBase_
