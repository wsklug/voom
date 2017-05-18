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
#include "HarmonicProteinBody.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{
  HarmonicProteinBody::HarmonicProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR, vector<vector<int > > HarmonicConn, double HarmonicStrength, double rHarmonic):
    ProteinBody(Proteins, Mat, SearchR),  _harmonicConn(HarmonicConn), _harmonicStrength(HarmonicStrength), _rHarmonic(rHarmonic)
  {
    
#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    
    // // Initialize Body.h containers
    // _nodes.insert(_nodes.begin(), Proteins.begin()->getHost(), Proteins.end()->getHost() );
    // for(ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
    //   _dof+=(*n)->dof();
    // }

  }; // ProteinBody constructor
	  

  
  //! Compute E0, E1, E2
  void HarmonicProteinBody::compute( bool f0, bool f1, bool f2 )
  {
    // Initialize energy to be zero
    if(f0) {
      _energy = 0.0;
    
      // Loop over protein elements
      for (uint i = 0; i < _prElements.size(); i++)
      {
	ProteinNode * A = _proteins[i];
	vector<ProteinNode *> domain = _prElements[i];
	for (uint j = 0; j < domain.size(); j++)
	{
	  _energy += _mat->computeEnergy(A, domain[j]);
	}
      }
      
      for (uint i = 0; i < _harmonicConn.size(); i++)
      {
	ProteinNode * A = _proteins[_harmonicConn[i][0]];
	ProteinNode * B = _proteins[_harmonicConn[i][1]];
	double r = A->getDistance(B);
	_energy += 0.5*_harmonicStrength*pow(r - _rHarmonic, 2.0);
      }
      
    } // if(f0) loop

    //To compute forces
    if(f1){
      // TODO
      // for (uint i=0; i < _prElements.size(); i++){
      // 	ProteinNode * A = _proteins[i];
      // 	vector<ProteinNode *> domain = _prElements[i];
      // 	for (uint j = 0; j < domain.size(); j++)
      // 	  {
      // 	    _mat->computeForce(A, domain[j]);
      // 	  }
      // }
    }
  };


} // namespace voom
