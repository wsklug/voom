// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------
//
// Reference:
//  Nonlinear quasi-Newton solver using L-BFGS code of Jorge Nocedal.
//  http://www.ece.northwestern.edu/~nocedal/lbfgs.html
//
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<cstdio>
#include<string>

#include "Lbfgs.h"

extern "C" void lbfgs_(int * n, int * m, 
		       double * x, double * f, double * g, 
		       int * diagco, double * diag,
		       int * iprint ,double * eps, 
		       double * xtol, double * w, 
		       int * iflag);

using namespace blitz;

using std::cout;
using std::endl;
using std::setw;
using std::right;
using std::left;
using std::scientific;

namespace voom
{
  void Lbfgs::resize(size_t n)  { 
    cout << "Lbfgs: resizing arrays for " << n << " DOF...";
    cout.flush();
    _x.resize(n);		
    _g.resize(n);
    _w.resize(2*_m*n+n+2*_m);
    _diag.resize(n);

    cout << " and initializing arrays to zero";
    cout.flush();
    _n = n;
    _f = 0.0;
    _x = 0.0;
    _g = 0.0;
    _w = 0.0;
    for(int idiag=0; idiag<n; idiag++) _diag(idiag) = 1.0;
    cout << '.' << endl;

  }      

  int Lbfgs::solve(Model * m) {

    _model = m;
    // resize arrays if necessary
    if( _n != _model->dof() || _x.size() != _model->dof() ) {
      resize( _model->dof() );
    }
    assert( _x.size() == _model->dof() );

    // copy starting guess from model
    _model->getField( *this );

    cout << "================================================================================"
	 << endl << endl
	 << "Starting (unbounded) BFGS iterations."
	 << endl << endl;

    cout << setw(14) << scientific << right << "|g|" 
	 << setw(14) << scientific  << right << "f" 
	 << setw(14) << right << "iterations" 
	 << endl
	 << "--------------------------------------------------------------------------------"
	 << endl;

    double xtol = 1.0e-16;

    _computeAll();
    cout << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	 << setw(14) << scientific  << right << _f
	 << setw(14) << right << 0  
	 << endl;

    int iflag = 0;

    int diagco = 0;

    lbfgs_(&_n,&_m,_x.data(),&_f,_g.data(),&diagco,_diag.data(),_print.data(),&_tol,&xtol,_w.data(),&iflag);

    int iters = 1;

    while(iflag==1) {
      _computeAll();
      if(_iprint > 0 && iters%_iprint==0) {
	cout << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	     << setw(14) << scientific  << right << _f
	     << setw(14) << right << iters  
	     << endl;
      }
      lbfgs_(&_n,&_m,_x.data(),&_f,_g.data(),&diagco,_diag.data(),_print.data(),&_tol,&xtol,_w.data(),&iflag);
      iters++;
      if(iters > _maxIterations) iflag = -10;
    }

    if(iflag==0) { // normal termination //
      _computeAll();
      cout << "LBFGS: Minimization program exited normally" << std::endl;
      cout << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	   << setw(14) << scientific  << right << _f
	   << setw(14) << right << iters  
	   << endl;
    }

    else { // bad termination, including iflag=-10 (maxIterations exceeded) //
      _computeAll();
      cout << "Something's gong wrong; iflag = " << iflag << std::endl;
      cout << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	   << setw(14) << scientific  << right << _f
	   << setw(14) << right << iters  
	   << endl;
    }

    cout.unsetf(ios_base::scientific);

    return 0;
  }



}; // end namespace
