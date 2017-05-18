// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__HarmonicProteinBody_h__)
#define __HarmonicProteinBody_h__

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

#include "ProteinBody.h"
#include "ProteinPotential.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{

  class HarmonicProteinBody : public ProteinBody
  {
  public:
    HarmonicProteinBody(vector<ProteinNode *> & Proteins, ProteinPotential * Mat, double SearchR, vector<vector<int > > HarmonicConn, double HarmonicStrength, double rHarmonic);
    
    //! Destructor
    ~HarmonicProteinBody() {};
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );

  private:

    // Springs positions
    vector<vector<int > > _harmonicConn;

    // Spring constant
    double _harmonicStrength;

    // Spring equilibrium distance
    double _rHarmonic;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  }; // Harmonic protein body

} // namespace voom

#endif // __HarmonicProteinBody_h__
