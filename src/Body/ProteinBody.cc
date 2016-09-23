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
  ProteinBody::ProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR, double Pressure):
    _proteins(Proteins), _mat(Mat), _searchR(SearchR), _pressure(Pressure)
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
      
      _prElements.push_back(domain);
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
      
      _prElements[i] = domain;
    }
    
  };
  
	  
  
  //! Compute E0, E1, E2
  void ProteinBody::compute( bool f0, bool f1, bool f2 )
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
      _energy -= _pressure*pow(_mat->getEquilibriumR(), 2.0);
    } // if(f0) loop

    //To compute forces
    if(f1){
      for (uint i=0; i < _prElements.size(); i++){
	ProteinNode * A = _proteins[i];
	vector<ProteinNode *> domain = _prElements[i];
	for (uint j = 0; j < domain.size(); j++)
	  {
	    _mat->computeForce(A, domain[j]);
	  }
      }
    }
  };



  double ProteinBody::computedWdEqPar()
  {
    double dWdEqPar = 0.0;
    // Loop over protein elements
    for (uint i = 0; i < _prElements.size(); i++)
    {
      ProteinNode * A = _proteins[i];
      vector<ProteinNode *> domain = _prElements[i];
      for (uint j = 0; j < domain.size(); j++)
      {
	dWdEqPar += _mat->computedWdEqPar(A, domain[j]);
      }
    }
    return dWdEqPar;
  };



  double ProteinBody::computeddWddEqPar()
  {
    double ddWddEqPar = 0.0;
    // Loop over protein elements
    for (uint i = 0; i < _prElements.size(); i++)
    {
      ProteinNode * A = _proteins[i];
      vector<ProteinNode *> domain = _prElements[i];
      for (uint j = 0; j < domain.size(); j++)
      {
	ddWddEqPar += _mat->computeddWddEqPar(A, domain[j]);
      }
    }
    return ddWddEqPar;
  }



  // reset Equilibrium parameter
  void ProteinBody::resetEquilibrium()
  {
    double a_prev = _mat->getEquilibriumParam(), a_new = 0.0;
    double a_initial = a_prev;
    double error = 1.0, tol = 1.0e-8;
    unsigned int iter = 0, MaxIter = 100;

    while (error > tol && iter < MaxIter)
    {
      a_new = a_prev - this->computedWdEqPar()/this->computeddWddEqPar();
      error = fabs(a_new-a_prev);
      iter++;
      a_prev = a_new;
      _mat->setEquilibriumParam(a_new);
    }
    
    cout << endl << "ProteinBody reset equilibrium" << endl;
    cout << "Iterations = " << iter << endl;
    cout << "Error      = " << error << endl;
    cout << "Previous EQ param. = " << a_initial << " New EQ param. = " << a_new << endl << endl;
  }

} // namespace voom
