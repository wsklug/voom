// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__ConcentratedForce_h__)
#define __ConcentratedForce_h__

#include "Node.h"
#include "Element.h"

namespace voom
{

  template<int N>
  class ConcentratedForce : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    ConcentratedForce(Node_t * node, const VectorND & f) 
      : _node(node), _force(f) { }

    void compute(bool f0, bool f1, bool f2) {
      if(f0) {
	_energy = - dot(_force,_node->point());
      }

      if(f1) {
	for(int i=0; i<N; i++) {
	  _node->addForce(i,-_force(i));
	}
      }

      return;
    }

    double force(int i) const {
      assert( i >= 0 && i < N );
      return ( -_force(i) );
    }
    
    //! set force vector
    void setForce(const VectorND & f) { _force = f; }

    //! set a component of force 
    void setForce(int i, double f) { 
      assert( i >= 0 && i < N );
      _force(i) = f; 
    }

  private:
    
    Node_t * _node;
    VectorND _force;
    
  };
};

#endif // __ConcentratedForce_h__
