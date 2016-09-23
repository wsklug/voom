// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Constraint_h__)
#define __Constraint_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"

namespace voom
{
  class Constraint
  {
  public:

    virtual void predict() = 0;
    virtual void correct() = 0;

    virtual double getForce(int a, int i) const {return 0;}

  };

  template<int N>
  class ConstrainNodes : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<N> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef typename DefNodeContainer::iterator DefNodeIterator;
    typedef typename DefNodeContainer::const_iterator ConstDefNodeIterator;   

      ConstrainNodes( const DefNodeContainer & nodes, 
		      const blitz::Array<int,1> & nbd,
		      const blitz::Array<double,1> & lo,
		      const blitz::Array<double,1> & hi) {

      _nbd.resize(nbd.size());
      _lo.resize(lo.size());
      _hi.resize(hi.size());
      _active.resize(nbd.size());
      _x.resize(nodes.size()*N);

      _nodes = nodes;
      _nbd = nbd;
      _hi = hi;
      _lo = lo;

    }

    virtual void predict() {

      for(int a=0; a<_nodes.size(); a++) {
	for(int b=0; b<N; b++){
	  _x(N*a+b) = _nodes[a]->getPoint(b);
	}
      }

      for(int i=0; i<_x.size(); i++){
	_active[i] = false;

	if (_nbd(i) == 1 && _x(i) <= _lo(i)){      //lower bond
	  _x(i) = _lo(i);
	  _active[i] = true;
	}
	else if (_nbd(i) == 3 && _x(i) >= _hi(i)){ //high bond
	  _x(i) = _hi(i);
	  _active[i] = true;
	}
	else if (_nbd(i) == 2){                  //both high and lower bonds
	  if (_x(i) <= _lo(i)){
	    _x(i) = _lo(i); _active[i] = true;
	  }
	  else if (_x(i) >= _hi(i)){
	    _x(i) = _hi(i); _active[i] = true;
	  }
	}
      }


      for(int a=0; a<_nodes.size(); a++) {
	for(int b=0; b<N; b++){
	  if(_active[N*a+b]) _nodes[a]->setPoint(b, _x(N*a+b));
	}
      }

      return;
    }

    virtual void correct() {

      for(int a=0; a<_nodes.size(); a++) {
	for(int b=0; b<N; b++){
	  if(_active[N*a+b]) _nodes[a]->setForce(b, 0.0);
	}
      }

      return;
    }



  private:

    DefNodeContainer _nodes;

    blitz::Array<int,   1> _nbd;
    blitz::Array<double,1> _lo;
    blitz::Array<double,1> _hi;
    blitz::Array<double,1> _x;

    std::vector<bool> _active;


  };

} // namespace voom

#endif // __Constraint_h__
