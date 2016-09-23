// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Ankush Aggarwal & William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__LMEtetBody_h__)
#define __LMEtetBody_h__

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
#include "LMEtet.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

// External LAPACK routine to invert matrices 
extern "C" void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
		       double* w, double* work, int* lwork, int* info );

namespace voom
{
  class LMEtetBody
  {
  public:
    
    //! Construct body from material, nodes, volumes, and LME parameters
    LMEtetBody(vector< DeformationNode<3>* > & defNodes, vector<vector<unsigned int > > & connT, 
	       double beta, double search_radius, unsigned int quad_order, double tol, unsigned int nItMax, 
	       double rho, double E, double nu);
    
    //! Destructor
    ~LMEtetBody() {};
    
    void computeNormalModes(double* K, double* eval, unsigned int ModesNumber, bool EgvFlag);
    
  private:
    vector< DeformationNode<3>* > & _defNodes;
    vector<vector<unsigned int > > & _connT;
    double _beta;
    double _search_radius;
    vector<vector<double > > _NatQP;
    vector<double > _weightQP;
    double _tol;
    unsigned int _nItMax;
    double _rho;
    double _Kmat[3][3][3][3];
  };  
}

#include "LMEtetBody.cc"

#endif // __LMEtetBody_h__
