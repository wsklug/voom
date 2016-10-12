// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.3  2005/08/22 22:28:28  klug
// Model::setField renamed putField
//
// Revision 1.2  2005/06/27 04:05:11  klug
// *** empty log message ***
//
// Revision 1.1  2005/05/23 18:05:35  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#if !defined(__SimulatedAnnealing_h__)
#define __SimulatedAnnealing_h__

#include "voom.h"
#include<iostream>
#include<iomanip>
#include<cstring>
#include<string>
#include<blitz/array.h>
#include<vector>
#include "Solver.h"

#if defined (_WIN32)
#if defined (_MSC_VER) || defined(__INTEL_COMPILER)
#define finite _finite
#endif
#endif

namespace voom
{

  /*!  A concrete class for nonlinear conjugate gradient solver for
    static equilibrium of a Finite Element model.
  */

  class SimulatedAnnealing : public Solver
  {
 
  public:
    
    enum TempSchedule { 
      LINEAR,           // T ~ alpha*t
      FAST,             // T ~ 1/t ("Fast Annealing")
      EXPONENTIAL       // T ~ alpha^t,  0 < alpha < 1
    };
    
    //! Default Constructor
    SimulatedAnnealing(bool debug=false) {
      _debug = debug;
      setParameters();
    }

    //! destructor
    virtual ~SimulatedAnnealing() {}

    //! overloading pure virtual function solve()
    int solve(Model * m);

    //
    // set method,  algorithm and parameters
    // for solving the nonlinear equations

    void setParameters
    (
     TempSchedule sched=FAST,
     const unsigned nSteps=1000, 
     const double T01=0,
     const double T02=0,
     const double finalTratio=1.0e-8,
     const unsigned printStride=100
     )
    {
      _schedule = sched;
      _nSteps = nSteps;
      _T01 = T01;
      _T02 = T02;
      _finalTratio = finalTratio;
      _printStride = printStride;
    }
    
    double & field(int i) {return _x(i);}
    double & function() {return _f;}
    double & gradient(int i) {
      std::cerr << "No stiffness in Simulate Annealing solver." << std::endl;
      exit(0);
    };    
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in Simulate Annealing solver." << std::endl;
      exit(0);
    };
    
    const double field(int i) const {return _x(i);}
    const double function() const {return _f;}
    const double gradient(int i) const {
      std::cerr << "No gradient in Simulate Annealing solver." << std::endl;
      exit(0);
    };
    const double hessian(int i, int j) const {
      std::cerr << "No stiffness in Simulate Annealing solver." << std::endl;
      exit(0);
    };

    double & hessian(int i) { return hessian(i,i);}
    const double hessian(int i) const { return hessian(i,i);}
    
    virtual double * field() {return _x.data();}
    virtual double * gradient() {return 0;}
        
    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _f=0.0;
    }
  
    void resize(size_t sz) { 
      _x.resize(sz); 
      _xSaved.resize(sz);
      _size = sz;
      _f = 0.0;
      _x = 0.0;
    }
      
  private:	

    typedef blitz::Array<double,1> _Vector;

    double _f;
    double _fSaved;
    _Vector _x;
    _Vector _xSaved;
    
    size_t _size;

    TempSchedule _schedule;

    double _T1;
    double _T2;
    double _T01;
    double _T02;
    double _finalTratio;

    unsigned _nSteps;
    
    unsigned _printStride;

    bool _debug;

    void _compute(Model * model) {
      model->putField( *this );
      model->computeAndAssemble(*this,true,false,false);
    };

    bool _changeState(Model * model);

    void _printState(Model * model, std::string name, _Vector & x) {
      cycleArrays(x,_x);
      model->putField( *this );
      model->print(name);
      cycleArrays(x,_x);
      model->putField( *this );      
      std::ofstream ofs("SAdata");
      ofs << _f << std::endl
	  << _fSaved << std::endl;
    }
  };
  
}; // namespace voom

#endif // __SimulatedAnnealing_h__
