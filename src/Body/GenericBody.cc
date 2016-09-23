// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file GenericBody.cc

  \brief GenericBody is a generic implemention for the Body interface,
  designed to call the compute functions of a generic list of Elements
  connected to a list of shared nodes.

*/

#include "GenericBody.h"

namespace voom
{

  //! Default Constructor
  GenericBody::GenericBody()
  {
    _output=paraview; 
    _energy=0.0;
  }

  GenericBody::GenericBody(ElementContainer & elems, NodeContainer & nodes) 
  {
    _elements = elems;
    _nodes = nodes;
    _dof = 0;
    for(int a=0; a<_nodes.size(); a++) _dof += _nodes[a]->dof();

    _output=paraview; 
    _energy=0.0;
  }
    
  //! Do mechanics on GenericBody
  void GenericBody::compute( bool f0, bool f1, bool f2 ) 
  {

    // Predictor/corrector approach for constraint
    for(ConstraintIterator c=_constraints.begin(); c!=_constraints.end(); c++) {
      (*c)->predict();
    }

    // Initialize energy and forces
    if(f0) _energy = 0.0;

    for(ElementIterator e=_elements.begin(); e!=_elements.end(); e++) {
      (*e)->compute( f0, f1, f2 );
     }

    if(f0) { // assemble total energy
      for(ElementIterator e=_elements.begin(); e!=_elements.end(); e++) {
	_energy += (*e)->energy();
      }
    }

    // correction for all constraints
    for(ConstraintIterator c=_constraints.begin(); c!=_constraints.end(); c++) {
      (*c)->correct();
    }
    
    return;
  }
  
} // namespace voom

