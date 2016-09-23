// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/05/25 02:09:36  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#include "DirectLinearSolver.h"

extern "C" void dgesv_( int * N, int * NRHS, double * A, int * LDA, 
			int *IPIV, double * B, int * LDB, int * INFO );

namespace voom {

  int DirectLinearSolver::solve(Model * m) {
    if(_size!=m->dof()) {
      resize(m->dof());
      zeroOutData(false,true,true);
    }
    m->computeAndAssemble(*this,false,true,true);
    // copy Df to x since routine will overwrite RHS with solution.
    _x = -_Df;

    int N = _x.size();
    int NRHS = 1;
    double * const A = _DDf.data();
    int LDA = N;
    int * const IPIV = _IPIV.data();
    double * const B = _x.data();
    int LDB = N;
    int INFO = 0;

    dgesv_( &N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO );

    if( INFO < 0 ) {
      std::cout << "dgesv Error: argument " << INFO << "had an illegal value."
		<< std::endl;
    } else if( INFO > 0 ) {
      std::cout << "dgesv Error: Matrix diagonal element" << INFO 
		<< "is zero.  Matrix is sigular."<< std::endl;
    } else {
      m->addField(*this);
    }
    
  }
}; // end namespace voom
