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
#include "ProteinBody.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{
  ProteinBody::ProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR):
    _proteins(Proteins), _mat(Mat), _searchR(SearchR) 
  {
    
#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    
    // // Initialize Body.h containers
    // _nodes.insert(_nodes.begin(), DefNodes.begin(), DefNodes.end() );
    // for(ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
    //   _dof+=(*n)->dof();
    // }
        
    // Initialize protein objects
    for(uint i = 0; i < Proteins.size(); i++)
    {
      ProteinNode * A = Proteins[i];
      vector<ProteinNode *> domain;
      // Find neighbors of A
      for(uint j = 0; j < Proteins.size(); j++)
      {
	if( i != j && A->getDistance(Proteins[j]) <= _searchR) {
	  domain.push_back(Proteins[j]);
	}
      }
      
      _elements.push_back(domain);
    }

  }; // ProteinBody constructor


  
  void ProteinBody::recomputeNeighbors(const double searchR) {
    _searchR = searchR;

    // Initialize protein objects
    for(uint i = 0; i < _proteins.size(); i++)
    {
      ProteinNode * A = _proteins[i];
      vector<ProteinNode *> domain;
      // Find neighbors of A
      for(uint j = 0; j < _proteins.size(); j++)
      {
	if( i != j && A->getDistance(_proteins[j]) <= _searchR) {
	  domain.push_back(_proteins[j]);
	}
      }
      
      _elements[i] = domain;
    }
    
  };
  
	  
  
  //! Compute E0, E1, E2
  void ProteinBody::compute( bool f0, bool f1, bool f2 )
  {
    // Initialize energy to be zero
    if(f0) {
      _energy = 0.0;
    
      // Loop over material objects
      for (uint i = 0; i < _elements.size(); i++)
      {
	ProteinNode * A = _proteins[i];
	vector<ProteinNode *> domain = _elements[i];
	for (uint j = 0; j < domain.size(); j++)
	{
	  _energy += _mat->computeEnergy(A, domain[j]);
	}
      }
    }
    
    return;
  };


} // namespace voom
