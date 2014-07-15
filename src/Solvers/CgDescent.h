// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------
//
// References:
//
// A new conjugate gradient method with guaranteed descent and an
// efficient line search, SIAM Journal on Optimization, 16 (2005),
// 170-192.
//
// CG_DESCENT, A conjugate gradient method with guaranteed descent,
// ACM Transactions on Mathematical Software, 32 (2006), 113-137.
//
// http://www.math.ufl.edu/~hager/
//
/////////////////////////////////////////////////////////////////////////

/*! 
  \file Lbfgsb.h

  \brief Interface to a concrete class for Conjugate Gradient solver for static
  equilibrium of a (nonlinear) Finite Element model.

*/

#if !defined(__CgDescent_h__)
#define __CgDescent_h__

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

  class CgDescent : public Solver
  {
 
  public:
    
    typedef blitz::Array<double,1> Vector_t;
    
    //! Default Constructor
    CgDescent(int n, double tol) : _tol(tol) { 
      resize(n);
    }
    //! destructor
    virtual ~CgDescent() {}
    
    double & field(int i) {return _x(i);}
    double & function() {return _f;}
    double & gradient(int i) {return _g(i);}
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in CgDescent solver." << std::endl;
      exit(0);
    }
    
   double field(int i) const {return _x(i);}
   double function() const {return _f;}
   double gradient(int i) const {return _g(i);}
   double hessian(int i, int j) const {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    };
    
    double & hessian(int i) { return hessian(i,i);}
    double hessian(int i) const { return hessian(i,i);}

    const blitz::Array<double,1> & gradient() const {return _g;}
    const blitz::Array<double,2> & hessian() const;

    double * field() {return _x.data();}
    double * gradient() {return _g.data();}

    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _f=0.0;
      if(f1) _g=0.0;
    }

    int size() const { return _x.size();}

    void resize(size_t n){ 
      _x.resize(n); 
      _g.resize(n); 
      _work.resize(4*n);
      _x = 0.0;
      _g = 0.0;
      _work = 0.0;
      _f = 0.0;
      return;
    }

    //! overloading pure virtual function solve()
    int solve(Model * m);

  private:	

    double _f;
    Vector_t _x;
    Vector_t _g;
    Vector_t _work;

    double _tol;

  };
  
}; // namespace voom

// Global model pointer needed to wrap C++ member functions for C code.
extern voom::Model * cgModel;
extern voom::CgDescent * cgSolver;

extern double myvalue( double *x );

extern void mygrad(double *g, double *x);


// The file cg_descent_c.parm must be present in the directory from
// which the code is executed.  It should have the following format.

// .1        delta        (Wolfe line search parameter)
// .9        sigma        (Wolfe line search parameter)
// 1.e-6     eps          (perturbation parameter for computing fpert)
// .66       gamma        (required decay factor in interval)
// 5.        rho          (interval growth factor used to get bracketing interval)
// .01       eta          (lower bound for cg's beta_k)
// .01       psi0         (factor used in starting guess for iteration 1)
// .1        psi1         (factor previous step multiplied by in QuadStep)
// 2.        psi2         (factor previous step is multipled by for startup)
// 1.e-12    QuadCutOff   (QuadStep if relative change in f > QuadCutOff)
// 0.e-12    StopFact     (factor multiplying starting |grad|_infty in StopRule)
// 1.e-3     AWolfeFac    (AWolfe = 0 => set AWolfe = 1 if |f-f0| < AWolfe_Fac*Ck)
// 1.        restart_fac  (restart cg in restart_fac*n iterations)
// 500.      maxit_fac    (terminate in maxit_fac*n iterations)
// 0.        feps         (stop when value change <= feps*|f|)
// .7        Qdecay       (used in Qk update: Qk = Qdecay*Qk + 1)
// 50        nexpand      (number of grow/shrink allowed in bracket)
// 50        nsecant      (number of secant steps allowed in line search)
// 1         PertRule     (0 => eps, 1 => eps*Ck)
// 1         QuadStep     (use initial quad interpolation in line search)
// 0         PrintLevel   0 (no print) 1 (intermediate results)
// 1         PrintFinal   0 (no print) 1 (print error messages, final error)
// 1         StopRule     1 (|grad|_infty <= max(tol,|grad|_0*StopFact) 0 (... <= tol*(1+|f|))
// 0         AWolfe       0 (Wolfe -- see AWolfeFac above) 1 (approx Wolfe)
// 0         Step         0 (no initial line search guess) 1 (guess in step arg)
// 0         debug        0 (no debugging) 1 (check for no increase in f)


#endif // __CgDescent_h__
