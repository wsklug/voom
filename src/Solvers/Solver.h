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
// Revision 1.2  2005/06/27 03:51:48  klug
// Issues with const.
//
// Revision 1.1  2005/05/25 02:15:29  klug
// Initial checkin.
//
//----------------------------------------------------------------------

/*! 
  \file Solver.h

  \brief Solver is the virtual base class for all Finite Element
  solvers, which implement the concept of evolving a Finite Element
  model, possibly achieving static equilibrium or advancing forward
  one time step in dynamics.

*/

#if !defined(__Solver_h__)
#define __Solver_h__

#include<blitz/array.h>
#include<vector>
#include "Model.h"

namespace voom
{

/*!    The virtual base class for all Finite Element
  solvers, which implement the concept of evolving a Finite Element
  model, possibly achieving static equilibrium or advancing forward
  one time step in dynamics.
*/

class Solver
{
 
public:

  //! Default Constructor
  Solver() {};
  
  virtual int solve(Model * m) = 0;

  virtual int size() const = 0;

  virtual double & field(int i) = 0;
  virtual double & function() = 0;
  virtual double & gradient(int i) = 0;
  virtual double & hessian(int i, int j) = 0;

  virtual double field(int i) const = 0;
  virtual double function() const = 0;
  virtual double gradient(int i) const = 0;
  virtual double hessian(int i, int j) const = 0;

  virtual double & hessian(int i) { return hessian(i,i);}
  virtual double hessian(int i) const { return hessian(i,i);}

  virtual void zeroOutData(bool f0, bool f1, bool f2) = 0;

  virtual void resize(size_t sz) = 0;

};

// struct for solver type storage
struct Storage : public Solver {

  double _E;
  blitz::Array<double,1> _x;
  blitz::Array<double,1> _DE;
  blitz::Array<double,2> _DDE;

  double & field(int i) {return _x(i);}
  double & function() {return _E;}
  double & gradient(int i) {return _DE(i);}
  double & hessian(int i, int j) {return _DDE(i,j);}
  double & hessian(int i) {return _DDE(i,i);}

  double field(int i) const {return _x(i);}
  double function() const {return _E;}
  double gradient(int i) const {return _DE(i);}
  double hessian(int i, int j) const {return _DDE(i,j);}
  double hessian(int i) const {return _DDE(i,i);}

  double * field() {return _x.data();}
  double * gradient() {return _DE.data();}

  void zeroOutData(bool f0, bool f1, bool f2) {
    if(f0) _E=0.0;
    if(f1) _DE=0.0;
    if(f2) _DDE=0.0;
  }
  
  void resize(size_t sz) { 
    _x.resize(sz); 
    _DE.resize(sz); 
    _DDE.resize(sz,sz);
    _x = 0.0;
    _DE = 0.0;
    _DDE = 0.0;
  }

  int size() const {return _x.size();}

  int solve(Model * m) {return 0;} 

};

}; // namespace voom
#endif // __Solver_h__
