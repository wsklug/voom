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
//  Nonlinear quasi-Newton solver using L-BFGS-B code of Jorge Nocedal.
//  http://www.ece.northwestern.edu/~nocedal/lbfgsb.html
//
/////////////////////////////////////////////////////////////////////////

#include<iostream>
#include<fstream>
#include<cstdio>
#include<string>

#include "Lbfgsb.h"

#if defined(_WIN32)
#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define setulb_ SETULB
#endif
#endif

extern "C" void setulb_(int * n, int *m, 
			double * x, double * l, double * u, 
			int * nbd, 
			double * f, double * g, 
			double * factr, double * pgtol,
			double * wa, 
			int * iwa,
			char * task, 
			int * iprint,  
			char * csave, 
			int * lsave, 
			int * isave, 
			double * dsave );

using namespace blitz;

using std::cout;
using std::endl;
using std::setw;
using std::right;
using std::left;
using std::scientific;

namespace voom
{
  void Lbfgsb::resize(size_t n)  { 
    cout << "Lbfgsb: resizing arrays for " << n << " DOF...";
    cout.flush();
    _x.resize(n);		
    _g.resize(n); 	
    _l.resize(n);
    _u.resize(n);
    _nbd.resize(n);
    _iwa.resize(3*n);
    _wa.resize(2*_m*n+4*n+12*_m*_m+12*_m);

    cout << " and initializing arrays to zero";
    cout.flush();
    _n = n;
    _f = 0.0;
    _x = 0.0;
    _g = 0.0;
    _l = 0.0;
    _u = 0.0;
    _nbd = 0;
    _iwa = 0;
    _wa = 0.0;
    cout << '.' << endl;

  }      

  void Lbfgsb::setBounds(const IntArray & nbd, 
			 const Vector_t & l, const Vector_t & u) {
    if(nbd.size() != _n || l.size() != _n || u.size() != _n ) {
      cout << "Lbfgsb::setBounds(): input arrays are incorrectly sized."
		<< endl;
      return;
    }

    _nbd = nbd;
    _l = l;
    _u = u;
    return;
  }


  int Lbfgsb::solve(Model * m) {

    _model = m;
    // resize arrays if necessary
    if( _n != _model->dof() || _x.size() != _model->dof() ) {
      resize( _model->dof() );
    }
    assert( _x.size() == _model->dof() );

    // copy starting guess from model
    _model->getField( *this );
      
    // set up misc. arrays and data
    char task[60], csave[60];
    for(int i=0; i<60; i++) task[i] = csave[i] = '\0';
      
    int lsave[4];
    for(int i=0; i<4; i++) lsave[i]=0;

    double dsave[29];
    for(int i=0; i<29; i++) dsave[i]=0.0;

    int isave[44];
    for(int i=0; i<44; i++) isave[i]=0;

    cout << "================================================================================"
	 << endl << endl
	 << "Starting BFGS iterations."
	 << endl << endl;

    cout << setw(14) << scientific << right << "|proj g|" 
	 << setw(14) << scientific << right << "|g|" 
	 << setw(14) << scientific  << right << "f" 
	 << setw(14) << right << "iterations" 
	 << setw(14) << right << "evaluations" 
	 << endl
	 << "--------------------------------------------------------------------------------"
	 << endl;

    // We start the iteration by initializing task.
 
    sprintf(task,"START");

    int iprint=-1;

    _computeAll();

    // ------- the beginning of the loop ----------
    while(true) {
      
      // This is the call to the L-BFGS-B code.
      setulb_(&_n, &_m, _x.data(), _l.data(), _u.data(), _nbd.data(), 
	      &_f, _g.data(), &_factr, &_pgtol, _wa.data(),_iwa.data(), 
	      &(task[0]), &iprint, &(csave[0]),
	      &(lsave[0]),&(isave[0]),&(dsave[0])); 
      if(strncmp(task,"FG",2)==0) {
	// The minimization routine has returned to request the
	// function f and gradient g values at the current x.

	_computeAll();

	// Go back to the minimization routine.
	continue;
      }

      else if( strncmp(task,"NEW_X",5)==0 ) {
	// stop if maximum number of iterations has been reached
	if ( _maxIterations > 0 && isave[29] > _maxIterations ) {
	  break;
	}

	// The minimization routine has returned with a new iterate,
	// and we have opted to continue the iteration.
	if ( _iprint>0 && isave[29]%_iprint == 0 ){
// 	  char name[100];sprintf(name,"%d",isave[29]);
	  _model->print("lbfgsbiter");
	  cout << setw(14) << scientific << right << dsave[12]
	       << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	       << setw(14) << scientific  << right << _f
	       << setw(14) << right << isave[29] 
	       << setw(14) << right << isave[33] 
	      << endl;

	}

	continue;
      }

      else if( strncmp(task,"CONV",4)==0 ) {
	_computeAll();
	//cout << task << endl;
	std::cout.write(task,49) << std::endl;
	//_model->print("lbfgsbconv");
	break;
      }

      else if( strncmp(task,"ABNORM",6)==0 ) {
	_computeAll();
	std::cout.write(task,49) << std::endl;
	//cout << task << endl;
	//_model->print("abnormal");
	break;
      }

      // If task is neither FG nor NEW_X we terminate execution.
      else {
	//cout << task << endl;
	std::cout.write(task,49) << std::endl;
	break;
      }

      // ---------- the end of the loop -------------
    }
    cout << setw(14) << scientific << right << dsave[12]
	 << setw(14) << scientific << right << blitz::max(blitz::abs(_g))
	 << setw(14) << scientific << right << _f
	 << setw(14) << right << isave[29] 
	 << setw(14) << right << isave[33] 
	 << endl << endl
	 << "================================================================================"
	 << endl << endl;

    cout.unsetf(ios_base::scientific);


    // ---------- the end of solve() -------------
    _projg = dsave[12];
    _iterNo = isave[29];

    _computeAll();

    return 0;
  }



}; // end namespace
