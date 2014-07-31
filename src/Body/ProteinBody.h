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
    ProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR);
    
    //! Destructor
    ~ProteinBody() {};
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );

    void recomputeNeighbors(double searchR);
    
    //! General printing of a Paraview file
    void printParaview(const string name) const {
      // Not implemented
    };

  private:
    vector<ProteinNode* > & _proteins;

    // Protein elements
    vector<vector<ProteinNode* > > _elements;
    
    // Protein material
    ProteinPotential * _mat;

    // Search radius
    double _searchR;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  }; // Protein body

} // namespace voom

#endif // __ProteinBody_h__
