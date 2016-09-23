// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__ProteinBody_h__)
#define __ProteinBody_h__

#include <blitz/array.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "Body.h"
#include "ProteinPotential.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{

  class ProteinBody : public Body
  {
  public:
    ProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR, double Pressure = 0.0);
    
    //! Destructor
    ~ProteinBody() {};
    
    //! Do mechanics on Body
    void virtual compute( bool f0, bool f1, bool f2 );
    double computedWdEqPar();
    double computeddWddEqPar();
    
    void recomputeNeighbors(double searchR);

    void resetEquilibrium();

    ProteinPotential * getPotential() { return _mat; };
    
    //! General printing of a Paraview file
    void printParaview(const string name) const {
      // Not implemented
    };

  protected:
    vector<ProteinNode* > & _proteins;

    // Protein elements
    vector<vector<ProteinNode* > > _prElements;
    
    // Protein material
    ProteinPotential * _mat;

    // Search radius
    double _searchR;

    // Pressure to be included when energy is computed
    double _pressure;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  }; // Protein body

} // namespace voom

#endif // __ProteinBody_h__
