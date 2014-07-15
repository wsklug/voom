// -*- C++ -*-
//----------------------------------------------------------------------
//
//                             Mo Bai
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------
#include "AffinityElement.h"

namespace voom
{
  //! Constructor
  AffinityElement::AffinityElement(const Shape_t & s, const NodeContainer & nodes ):_s(s)
  {
    // _s = s;
    //! initialize NodeContainer
    unsigned nNodes = nodes.size();
      
    _nodes = nodes;
    
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
      _baseNodes.push_back(*n);
  } 

  const Vector2D & AffinityElement::Centroid() {
    return _centroid;
  }
    
  // return strain tensor
  const Tensor2D & AffinityElement::Strain() const {
    return _strain;
  }
  
  // return displacement gradient tensor
  const Tensor2D & AffinityElement::DisplacementGradient() const {
    return _dispGrad;
  }
  
  double AffinityElement::Rotation() const {
    return _rotation;
  }

  void AffinityElement::compute() {
      
    // Compute spatial derivatives of shape functions from
    // parametric derivatives by transforming with matrix jacobian
    // of isoparametric mapping.
    
    // parametric derivatives from shape object
    const Shape_t::DerivativeContainer & dnds = _s.derivatives();
      
    // matrix duds;
    tvmet::Matrix<double,2,2> duds(0.0);
    for(int i=0; i<2; i++) {
      for(int alpha=0; alpha<2; alpha++) {   
	for(int a=0; a<_nodes.size(); a++) {
	  duds(i,alpha) += dnds[a](alpha)*(_nodes[a]->getPoint(i) - _nodes[a]->getPosition(i) );
	}
      }      
    }
      
    // matrix jacobian dxds;
    tvmet::Matrix<double,2,2> dxds(0.0);
    for(int i=0; i<2; i++) {
      for(int alpha=0; alpha<2; alpha++) {   
	for(int a=0; a<_nodes.size(); a++) {
	  dxds(i,alpha) += dnds[a](alpha)*(_nodes[a]->getPosition(i) );
	}
      }      
    }
    
    // compute scalar jacobian by calculating the determinant of dxds
    double J = dxds(0,0)*dxds(1,1) - dxds(0,1)*dxds(1,0);
      
    // invert matrix jacobian
    tvmet::Matrix<double,2,2> invJac(0.0);
    invJac(0,0) = dxds(1,1);
    invJac(0,1) = -dxds(0,1);
    invJac(1,0) = -dxds(1,0);
    invJac(1,1) = dxds(0,0);
    
    invJac /= J;
    
    _dispGrad = duds * invJac;
    _strain = 0.5*(tvmet::trans(_dispGrad) + _dispGrad);
    _rotation = (_dispGrad(0,1)-_dispGrad(1,0))/2.0;

    //calculate the area of element
    assert(_nodes.size()==3);
    tvmet::Vector<double,2> AB, AC;
    AB = _nodes[1]->position() - _nodes[0]->position();
    AC = _nodes[2]->position() - _nodes[0]->position();
    _area = (AB(0)*AC(1) - AB(1)*AC(0))*0.5; 
    _centroid = 0.0,0.0;
    for(int k=0; k<3; k++) _centroid += _nodes[k]->position();
    _centroid /= 3.0;
  }
    
  // Access the container of nodes
  //const NodeContainer & AffinityElement::nodes() { return _nodes; }

} // namespace voom
