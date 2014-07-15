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
// Revision 1.9  2005/08/20 01:28:32  klug
// acinclude.m4
//
// Revision 1.8  2005/06/27 03:52:53  klug
// Issues with const.
//
// Revision 1.7  2005/05/23 18:03:17  klug
// Enforced consistency with new Solver interface.  Minor tweaks.
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
  \file ConjugateGradientWSK.h

  \brief Interface to a concrete class for nonlinear conjugate
  gradient solver for static equilibrium of a Finite Element model.

*/

#if !defined(__ConjugateGradientWSK_h__)
#define __ConjugateGradientWSK_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<vector>
#include "Solver.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*!  A concrete class for nonlinear conjugate gradient solver for
    static equilibrium of a Finite Element model.
  */

  class ConjugateGradientWSK : public Solver
  {
 
  public:
    
    // Line Search enum
    enum LineSearchChoice {
      Newton,
      Secant,
      Wolfe
    };
	
    // Algorithm enum
    enum AlgorithmChoice {
      SD,
      PR,
      FR
    };
    
    //! Default Constructor
    ConjugateGradientWSK(bool debug=false) {
      _debug = debug;
      setParameters();
      setWolfeParameters();
#ifdef WITH_MPI
      MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
      MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    }

    //! destructor
    virtual ~ConjugateGradientWSK() {}

    //! overloading pure virtual function solve()
    int solve(Model * m);

    //
    // set method,  algorithm and parameters
    // for solving the nonlinear equations

    void setParameters
    (
     LineSearchChoice lsMethod = Newton,			
     AlgorithmChoice  cgMethod = PR,
     int maxIter = 1000,
     int restartStride = 10,
     int printStride = 10,
     double tol = 1.0e-6,
     double absTol = 1.0e-6,
     double tolLS = 1.0e-3,
     int maxIterLS = 20
     )
    {
      _lsMethod = lsMethod;
      _cgMethod = cgMethod;
      _maxIter = maxIter;
      _restartStride = restartStride;
      _printStride = printStride;
      _tol = tol;
      _absTol = absTol;
      _tolLS = tolLS;
      _maxIterLS = maxIterLS; 
    }
    
    void setWolfeParameters
    (
     double  c1 = 1.0e-4,
     double  c2 = 0.9,
     double  min_alphaInc = 1.0
     )
    {
      _c1 = c1;
      _c2 = c2;
      _min_alpha_Inc = min_alphaInc;
    }
    
    double & field(int i) {return _x(i);}
    double & function() {return _f;}
    double & gradient(int i) {return _gradf(i);}
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    }
    
    double field(int i) const {return _x(i);}
    double function() const {return _f;}
    double gradient(int i) const {return _gradf(i);}
    double hessian(int i, int j) const {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    };
    
    double & hessian(int i) { return hessian(i,i);}
    double hessian(int i) const { return hessian(i,i);}

    const blitz::Array<double,1> & gradient() const {return _gradf;}
    const blitz::Array<double,2> & hessian() const;

    double * field() {return _x.data();}
    double * gradient() {return _gradf.data();}

    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _f=0.0;
      if(f1) _gradf=0.0;
    }
  
    void resize(size_t sz) { 
      _x.resize(sz);		
      _gradf.resize(sz); 	
      _size = sz;
      _f = 0.0;
      _x = 0.0;
      _gradf = 0.0;
      assert(_x.size()==sz); 	  
      assert(_gradf.size()==sz);
    }
      
    int size() const { return _size;}

  private:	

    typedef blitz::Array<double,1> Vector_t;

    double _f;
    Vector_t _x;
    Vector_t _gradf;
    
    size_t _size;

    Model * _model;

    // using which line search method
    LineSearchChoice _lsMethod;

    //  using PR or FR algorithm
    AlgorithmChoice _cgMethod;

    int _maxIter;

    bool _debug;

    double _tol;

    double _tolLS;

    double _absTol;

    int _printStride;

    int _restartStride;

    // additional parameters for Wolfe line search method
    // see reference for meaning
    double _c1;
    double _c2;
    double _min_alpha_Inc;
    int    _maxIterLS;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif

    double _computeAll( Vector_t & x, Vector_t & grad );
//     //! update f
    double _computeFunction( Vector_t & x ) ;
//     double _computeFunction( const Vector_t & x ) {
//       for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
// 	if( fabs( *i ) > 10.0 ) {
// 	  std::cerr << *i << std::endl
// 		    << x << std::endl;
// 	  exit(0);
// 	}
//       }
//       _model->setPositions( x );
//       _model->compute( true, false, false );
//       return _model->getEnergy();
//     }

//     //! update gradf
    void _computeGradient( Vector_t & x, Vector_t & grad );
//     void _computeGradient( const Vector_t & x, Vector_t & grad ) {
//       for( Vector_t::const_iterator i=x.begin(); i!=x.end(); ++i ) {
// 	if( fabs( *i ) > 10.0 ) {
// 	  std::cerr << *i << std::endl
// 		    << x << std::endl;
// 	  exit(0);
// 	}
//       }
//       _model->setPositions( x );
//       _model->compute( false, true, false );
//       _model->getResidual( grad );
//       grad = -grad;
//       return;
//     }

    double _computeDirectionalD(Vector_t& x, 
				const Vector_t & dir,
				Vector_t & grad  ) {
      _computeGradient( x, grad );
      return _directionalD( dir, grad );
    };

    inline 
    double _directionalD(const Vector_t & dir, const Vector_t & grad) const {
      return blitz::sum( grad*dir );
    }

    //
    // Line search algorithm satisfying strong Wolfe conditions
    //
    double _lineSearch(const Vector_t & dir);
    double _secantLineSearch(const Vector_t & dir);
    double _wolfeLineSearch(const Vector_t & dir);
    
    //
    // Specialty function used by Wolfe line search algorithm
    //
    double _zoom(const double c1, const double c2, 
		 const double f0, const double df0, 
		 const Vector_t & dir, 
		 Vector_t & x, Vector_t & grad,
		 double lo, double hi);

  };
  
}; // namespace voom

#endif // __ConjugateGradientWSK_h__
