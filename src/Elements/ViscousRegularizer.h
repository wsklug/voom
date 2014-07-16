// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__ViscousRegularizer_h__)
#define __ViscousRegularizer_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace voom {

  //! Simple viscous regularization
  /*! ViscousRegularizer computes an energy quadratic in the difference
    between the current vector of DOF of a set of nodes and those at
    some previous state/increment/iteration (termed the reference
    state).  
    This difference can be thought of as a "velocity".  The
    force is thus proportional to this velocit, hence the label
    "viscous".

    \f[
    E = \frac{k}{2}\sum_{ia} (x_{ia} - \bar{x}_{ia})^2
    \f]

    \f[
    f_{ia} = k(x_{ia}-\bar{x}_{ia})
    \f]

    \f[
    k_{iajb} = k\delta_{ab}\delta_{ij}
    \f]
  */
  class ViscousRegularizer :
    public Element {

    public:
    using Element::BaseNodeContainer; 
    using Element::BaseNodeIterator; 
    using Element::ConstBaseNodeIterator; 

    //! Construct from a set of nodes and a "viscosity" proportionality factor
    ViscousRegularizer( const BaseNodeContainer & nodes, double viscosity ) 
      {
	_baseNodes = nodes;
	_viscosity = viscosity;
	_energy = 0.0;
	int nNodes = _baseNodes.size();
	int maxdof = 0;
	int dof=0;
// 	for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
// 	  dof+=(*n)->dof();
// 	}
	for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
	  if((*n)->dof() > maxdof) maxdof = (*n)->dof();
	}
	//	_reference.resize(dof);
	_reference.resize(nNodes,maxdof);
	_reference = 0.0;
	step();
      }

    //! Assign the current state to the reference state
    void step() {
//       int I = 0;
//       for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
// 	for(int i=0; i<(*n)->dof(); i++,I++) {
// 	  _reference(I) = (*n)->getPoint(i);
// 	}
//       }

      int nNodes = _baseNodes.size();

      // parallel block //
#ifdef _OPENMP	
#pragma omp parallel default(shared)
#endif
      {
#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif
	for(int n=0; n<nNodes; n++) {
	  NodeBase* node = _baseNodes[n];
	  int ndof = node->dof();
	  for(int i=0; i<ndof; i++) {
	    _reference(n,i) = node->getPoint(i);
	  }
	}
      }
      // end parallel block
    }

    //! Compute the energy, force, and stiffness
    void compute(bool f0, bool f1, bool f2) {
      if(f0) _energy = 0.0;

//       int I=0;
//       for(BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
// 	for(int i=0; i<(*n)->dof(); i++,I++) {
// 	  double dx = (*n)->getPoint(i) - _reference(I);
// 	  if(f0) _energy += 0.5 * _viscosity * dx * dx;
// 	  if(f1) (*n)->addForce(i, _viscosity * dx); 
// 	}
//       }

      // parallel block //
      int nNodes = _baseNodes.size();
      double tmpenergy = 0.0;
#ifdef _OPENMP	
#pragma omp parallel default(shared)
#endif
      {
#ifdef _OPENMP
#pragma omp for schedule(static) nowait reduction(+:tmpenergy)
#endif
	for(int n=0; n<nNodes; n++){
	  NodeBase* node = _baseNodes[n];
	  int ndof = node->dof();
	  for(int i=0; i<ndof; i++) {
	    double dx = node->getPoint(i) - _reference(n,i);
	    if(f0) tmpenergy += 0.5 * _viscosity * dx * dx;
	    if(f1) node->addForce(i, _viscosity * dx); 
	  }
	}
      }
      // end parallel block //

      _energy = tmpenergy;
      
    }

    //! return the viscosity proportionality factor
    double viscosity() const {return _viscosity;}

    //! return the viscosity proportionality factor
    void setViscosity(double v) {_viscosity = v;}

    //! compute norm of difference between reference and current
    double velocity() const {
      double v=0.0;
      int I=0;
      for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
	for(int i=0; i<(*n)->dof(); i++,I++) {
	  v = std::max( v, std::abs( _reference(I) - (*n)->getPoint(i) ) );
	}
      }
      return v;
    }
    private:

    //! proportionality factor \f$k\f$ in the energy
    double _viscosity;

    //! reference state \f$\bar{x}_{ia}\f$
    blitz::Array<double, 2> _reference;
  };

} // end namespace

#endif // __ViscousRegularizer_h__
