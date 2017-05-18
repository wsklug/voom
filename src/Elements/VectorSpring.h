// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__VectorSpring_h__)
#define __VectorSpring_h__

#include "Node.h"
#include "Element.h"

namespace voom
{

  template<int N>
  class VectorSpring : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    VectorSpring(Node_t * node,  double k) 
      : _node(node), _k(k) { target = false; }

    void compute(bool f0, bool f1, bool f2) {
      const VectorND & X = _node->position();
      const VectorND & x = _node->point();
      
      

      double v0;
      if (target) 
	v0 = _targetv0;
      else
	v0 = norm2(X);

      double v  = norm2(x);

      if(f0) {
	_energy = 0.5*_k*sqr(v-v0);
      }

      if(f1) {	
	for(int i=0; i<N; i++) {	  
	  double f = _k*(v-v0)*x(i)/v;
	  _node->addForce(i, f);	  
	}
      }

      return;
    }
    
    double stiffness() const {return _k;}

    void setStiffness(double k) { _k = k; }

    void setTargetv0 (double v0) {_targetv0 = v0;}

    void setTarget(bool flag) {target = flag;}

    double vectorMag() const {
      const VectorND & x = _node->point();
      
      double v  = norm2(x);
      return v;
    }

  private:
    
    Node_t * _node;
    double _k;
    double _targetv0;
    bool target;
    
  };
};

#endif // __Spring_h__
