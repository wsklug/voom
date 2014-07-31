// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__ViscosityBody_h__)
#define __ViscosityBody_h__

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
#include "LennardJones.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{
  class ViscosityBody : public Body
  {
  public:
    ViscosityBody(const vector<DeformationNode<3> *> & DefNodes,
		  const vector<tvmet::Vector<int, 3> > & Connectivity,
		  double sigma);
    
    //! Destructor
    ~ViscosityBody() {};
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );

    void resetRefConf();

    //! Return the energy of the body
    double totalStrainEnergy() const { return _energy; };
    
    //! General printing of a Paraview file
    void printParaview(const string name) const {
      // Not implemented
    };

    void setSigma(double NewSigma) { _sigma = NewSigma; };
    


  private:
    // List of 3D nodes pairs
    vector<vector<DeformationNode<3> *> > SpringsNodes;

    // List of spring lenghts
    vector<double > SpringsLengths;

    // Nodes
    const vector<DeformationNode<3> *> & _defNodes;

    // Spring constant
    double _sigma;
    
#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  }; // Viscosity body

} // namespace voom

#endif // __ViscosityBody_h__
