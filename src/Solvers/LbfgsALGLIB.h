// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------
//
/////////////////////////////////////////////////////////////////////////

/*! 
  \file LbfgsALGLIB.h

  \brief Interface to a concrete class for L-BFGS solver for static
  equilibrium of a (nonlinear) Finite Element model.

*/

#if !defined(__LbfgsALGLIB_h__)
#define __LbfgsALGLIB_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<vector>
#include "Solver.h"
#include "lbfgs.h"


namespace voom
{

  /*!  A concrete class for nonlinear conjugate gradient solver for
    static equilibrium of a Finite Element model.
  */

  class LbfgsALGLIB : public Solver
  {
 
  public:
    
    typedef blitz::Array<double,1> Vector_t;
    typedef blitz::Array<int,1> IntArray;
    
    //! Default Constructor
    LbfgsALGLIB(int m=5, 
	   double EpsF=1.0e-6, double EpsG=1.0e-6, double EpsX=1.0e-10, 
	   int iprint=0, int maxIterations=0) 
      : _iprint(iprint) 
    {
      _lbstate.m = m;
      _lbstate.epsg = EpsG;
      _lbstate.epsf = EpsF;
      _lbstate.epsx = EpsX;
      _lbstate.maxits = maxIterations;
      _lbstate.flags = 0;
    }

    //! destructor
    virtual ~Lbfgsb() {}
    
    double & field(int i) {return _lbstate.x(i);}
    double & function() {return _lbstate.f;}
    double & gradient(int i) {return _lbstate.g(i);}
    
    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _lbstate.f=0.0;
      if(f1) _lbstate.g=0.0;
    }

    int size() const { return _n;}

    void resize(size_t n);

    //! overloading pure virtual function solve()
    int solve(Model * m);

    int iterationNo();

  private:	
 
    lbfgsstate _lbstate;
    lbfgsreport _lbreport;

    int _iprint;

    Model * _model;

    void _computeAll() {
      _model->putField( *this );
      _model->computeAndAssemble( *this, true, true, false );
      return;
    }
    
  };
  
}; // namespace voom

#endif // __LbfgsALGLIB_h__

