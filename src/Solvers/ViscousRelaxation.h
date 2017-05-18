// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2005 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__ViscousRelaxation_h__)
#define __ViscousRelaxation_h__

#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<vector>
#include "Solver.h"

namespace voom
{

  /*!  A concrete class for a viscous relaxation solver for
    static equilibrium of a Finite Element model.
  */

  class ViscousRelaxation : public Solver
  {
 
  public:

    typedef blitz::Array<double,1> Vector_t;
    
    ViscousRelaxation(int n,
		      double dt=1.0e-8, 
		      double tol=1.0e-6,
		      double absTol=1.0e-6,
		      int maxIter=1000,
		      int printStride=100,
		      bool debug=false) 
      : _dt(dt), _tol(tol), _absTol(absTol), 
	_maxIter(maxIter), _printStride(printStride), _debug(debug) 
    {
      resize(n);
    }

    //! destructor
    virtual ~ViscousRelaxation() {}

    //! overloading pure virtual function solve()
    int solve(Model * m);

    double & field(int i) {return _x(i);}
    double & function() {return _f;}
    double & gradient(int i) {return _gradf(i);}
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      exit(0);
    }
    
    const double field(int i) const {return _x(i);}
    const double function() const {return _f;}
    const double gradient(int i) const {return _gradf(i);}
    const double hessian(int i, int j) const {
      std::cerr << "No stiffness in ConjugateGradient solver." << std::endl;
      return 0;
    };

    double & hessian(int i) { return hessian(i,i);}
    const double hessian(int i) const { return hessian(i,i);}
    
//     const blitz::Array<double,1> & gradient() const {return _gradf;}
//     const blitz::Array<double,2> & hessian() const;

    double * field() { return _x.data();}
    double * gradient() { return _gradf.data();}
    
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

    Vector_t _x;
    Vector_t _gradf;

    double _f;
    double _dt;
    double _tol;
    double _absTol;
    
    size_t _size;

    int _maxIter;
    int _printStride;

    bool _debug;

    Model * _model;

    //! update gradf
    void _computeGradient( Vector_t & x, Vector_t & grad );

  };
  
}; // namespace voom

#endif // __ViscousRelaxation_h__
