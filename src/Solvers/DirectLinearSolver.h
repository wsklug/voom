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
// Revision 1.2  2005/06/27 03:52:53  klug
// Issues with const.
//
// Revision 1.1  2005/05/25 02:09:36  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#if !defined(__DirectLinearSolver_h__)
#define __DirectLinearSolver_h__

#include<blitz/array.h>
#include<vector>
#include "Model.h"
#include "Solver.h"

namespace voom
{

/*! Wrapper class using Lapack routines for direct solution of linear
  equations of the form Ax=b.
*/

  class DirectLinearSolver : public Solver
{
 
public:

  //! Default Constructor
  DirectLinearSolver() {};
  
  int solve(Model * m);

  double & field(int i) {return _x(i);}
  double & function() {return _f;}
  double & gradient(int i) {return _Df(i);}
  double & hessian(int i, int j) {return _DDf(i,j);}
  double & hessian(int i) {return _DDf(i,i);}

  double const field(int i) const {return _x(i);}
  double const function() const {return _f;}
  double const gradient(int i) const {return _Df(i);}
  double const hessian(int i, int j) const {return _DDf(i,j);}
  double const hessian(int i) const {return _DDf(i,i);}

//   const blitz::Array<double,1> & gradient() const {return _Df;}
//   const blitz::Array<double,2> & hessian() const {return _DDf;}
  double * field() { return _x.data();}
  double * gradient() { return _Df.data();}

  void zeroOutData(bool f0, bool f1, bool f2) {
    if(f0) _f=0.0;
    if(f1) _Df=0.0;
    if(f2) _DDf=0.0;
  }
  
  void resize(size_t sz) { 
    _x.resize(sz); 
    _Df.resize(sz); 
    _DDf.resize(sz,sz);
    _IPIV.resize(sz);
    _size = sz;
}

 private:
  size_t _size;

  double _f;

  blitz::Array<double,1> _x;
  blitz::Array<double,1> _Df;
  blitz::Array<double,2> _DDf;

  blitz::Array<int,1> _IPIV;
};

}; // namespace voom
#endif // __DirectLinearSolver_h__
