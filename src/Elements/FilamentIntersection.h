// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__IntersectionFinder_h__)
#define __IntersectionFinder_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"

namespace voom
{

  template<int N>
  class IntersectionFinder : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    Spring(Node_t * nodeA, Node_t * nodeB, double k) 
      : _nodeA(nodeA), _nodeB(nodeB), _k(k) { 
      const VectorND & XA = _nodeA->position();
      const VectorND & XB = _nodeB->position();
      _d0 = norm2(XB-XA);

      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);
    }

    Spring(Node_t * nodeA, Node_t * nodeB, double k, double d) 
      : _nodeA(nodeA), _nodeB(nodeB), _k(k), _d0(d) {
      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);
    }

    void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double d  = norm2(xB-xA);

      if(f0) {
	_energy = 0.5*_k*sqr(d-_d0);
      }

      if(f1) {	
	for(int i=0; i<N; i++) {	  
	  double f = _k*(d-_d0)*(xA(i)-xB(i))/d;
	  _nodeA->addForce(i, f);
	  _nodeB->addForce(i,-f);
	}
      }

      return;
    }
    
    double stiffness() const {return _k;}

    void setStiffness(double k) { _k = k; }

    void resetLength() {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
       _d0  = norm2(xB-xA);
       return;
    }
  private:
    
    Node_t * _nodeA;
    Node_t * _nodeB;
    double _k;
    double _d0;
    
  };
};

#endif // __Spring_h__
