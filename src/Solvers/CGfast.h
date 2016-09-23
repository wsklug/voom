// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   Feng Feng and William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/08/20 01:28:32  klug
// acinclude.m4
//
// Revision 1.1  2005/05/23 18:07:36  klug
// Initial checkin.
//
//----------------------------------------------------------------------
//
//    Reference:
//
//       Jonathan Richard Shewchuk, "An introduction to the conjugate
//       gradient method without agonizing pain"
//       Edition 5/4:  August 4, 1994
//       School of Computer Science, Carnegie Mellon University
//       Pittsburgh, PA 15213
//
//  http://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.ps
//
/////////////////////////////////////////////////////////////////////////

/*! 
  \file CGfast.h

  \brief Interface to a concrete class for nonlinear conjugate
  gradient solver for static equilibrium of a Finite Element model.

*/

#if !defined(__CGfast_h__)
#define __CGfast_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<vector>
#include "Solver.h"

namespace voom
{

  /*!  A concrete class for nonlinear conjugate gradient solver for
    static equilibrium of a Finite Element model.
  */

  class CGfast : public Solver
  {
 
  public:
    
    // Algorithm enum
    enum AlgorithmChoice {
      SD,
      PR,
      FR
    };
    
    //! Default Constructor
    CGfast(bool debug=false) {
      _debug = debug;
      setParameters();
    }

    //! destructor
    virtual ~CGfast() {}

    //! overloading pure virtual function solve()
    int solve(Model * m);

    //
    // set method,  algorithm and parameters
    // for solving the nonlinear equations

    void setParameters
    (
     AlgorithmChoice  cgMethod = PR,
     int maxIter = 1000,
     int restartStride = 10,
     int printStride = 10,
     double tol = 1.0e-6,
     double absTol = 1.0e-6,
     double tolLS = 1.0e-3,
     int maxIterLS = 20,
     double sigma = 1.0e-6
     )
    {
      _cgMethod = cgMethod;
      _maxIter = maxIter;
      _restartStride = restartStride;
      _printStride = printStride;
      _tol = tol;
      _absTol = absTol;
      _tolLS = tolLS;
      _maxIterLS = maxIterLS; 
      _sigma = sigma;
    }
    
    double & field(int i) {return _x(i);}
    double & function() {return _f;}
    double & gradient(int i) {return _g(i);}
    double & hessian(int i) {return _h(i);}
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    }
    
    const double field(int i) const {return _x(i);}
    const double function() const {return _f;}
    const double gradient(int i) const {return _g(i);}
    const double hessian(int i) const { return _h(i);}
    const double hessian(int i, int j) const {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    };
    
    const blitz::Array<double,1> & gradient() const {return _g;}
    const blitz::Array<double,2> & hessian() const;
    
    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _f=0.0;
      if(f1) _g=0.0;
      if(f2) _h=0.0;
    }
  
    void resize(size_t sz) { 
      _x.resize(sz); 
      _g.resize(sz); 
      _h.resize(sz); 
      _size = sz;
      _f = 0.0;
      _x = 0.0;
      _g = 0.0;
      _h = 0.0;
    }

    int size() const {return _size;}

  private:	

    typedef blitz::Array<double,1> Vector_t;
    
    //! function value
    double _f;   

    //! argument
    Vector_t _x;

    //! gradient
    Vector_t _g;

    //! diagonal of Hessian
    Vector_t _h;
    
    size_t _size;

    Model * _model;

    //  using PR or FR algorithm
    AlgorithmChoice _cgMethod;

    int _maxIter;

    bool _debug;

    double _tol;

    double _tolLS;

    double _absTol;

    double _sigma;

    int _printStride;

    int _restartStride;

    int    _maxIterLS;

    void _compute(bool f0, bool f1, bool f2=false) ;

  };
  
}; // namespace voom

#endif // __CGfast_h__
