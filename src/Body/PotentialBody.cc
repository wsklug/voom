// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
#include <string>
#include <fstream>
#include <blitz/array-impl.h>
#include "VoomMath.h"
#include "PotentialBody.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{
  PotentialBody::PotentialBody(Potential * Mat,const vector<DeformationNode<3> * > & DefNodes,
                 	       double SearchR):
	_mat(Mat), _defNodes(DefNodes), _searchR(SearchR) 
  {
    
#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    
    // Initialize Body.h containers
    _nodes.insert(_nodes.begin(), DefNodes.begin(), DefNodes.end() );
    for(ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
      _dof+=(*n)->dof();
    }
        
    // Initialize material objects
    _elementVector.reserve(_defNodes.size());

    for(uint i =0; i < _defNodes.size(); i++)
    {
      Vector3D CenterNode = _defNodes[i]->point();
      set<DeformationNode<3> *> domain;
      // Find neighbors of CenterNode
      for(uint j = 0; j < _defNodes.size(); j++)
      {
	if( tvmet::norm2(CenterNode - _defNodes[j]->point()) <= _searchR && i != j) {
	  domain.insert(_defNodes[j]);
	}
      }
      // Build potential element
      PotentialElement * el = new PotentialElement(_mat, _defNodes[i], domain);
      _elementVector.push_back(el);
    }

  }; // PotentialBody constructor


  
  void PotentialBody::recomputeNeighbors(const double searchR) {
    _searchR = searchR;
    
    for(uint i =0; i < _defNodes.size(); i++)
    {
      Vector3D CenterNode = _defNodes[i]->point();
      set<DeformationNode<3> *> domain;
      // Find neighbors of CenterNode
      for(uint j = 0; j < _defNodes.size(); j++)
      {
	if( tvmet::norm2(CenterNode - _defNodes[j]->point()) <= _searchR && i != j) {
	  domain.insert(_defNodes[j]);
	}
      }
      // Modify domain of each element
      _elementVector[i]->resetDomain(domain);

    }
    
  };
  
	  
  
  //! Compute E0, E1, E2
  void PotentialBody::compute( bool f0, bool f1, bool f2 )
  {
    // Initialize energy to be zero
    if(f0) _energy = 0.0;
    
    // Loop over material objects
    for (uint i = 0; i < _elementVector.size(); i++)
    {
      _elementVector[i]->compute(f0, f1, f2);
      _energy += _elementVector[i]->energy();
    }

    return;
  };


} // namespace voom
