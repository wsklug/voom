// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file FESolver.h

  \brief FESolver is the virtual base class for all Finite Element
  solvers, which implement the concept of evolving a Finite Element
  model, possibly achieving static equilibrium or advancing forward
  one time step in dynamics.

*/

#if !defined(__FESolver_h__)
#define __FESolver_h__

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
//template < class key_t >
class FESolver
{
 
public:

//   typedef std::pair< Node*, int > DOF;
  
  //! Default Constructor
  FESolver() {};
  
  virtual int solve() = 0;

  virtual void setModel( Model * const m ) = 0;

//   // Pass field values between solver and its client; const functions
//   // used on right of the = operator (getting values from solver),
//   // non-const functions on the left (setting values in solver).

//   //! positions
//   const double 	positions( const DOF & dof ) const = 0;
//   double & 	positions( const DOF & dof ) = 0;

//   //! positions
//   const double 	velocities( const DOF & dof ) const = 0;
//   double & 	velocities( const DOF & dof ) = 0;

//   //! positions
//   const double 	accelerations( const DOF & dof ) const = 0;
//   double & 	accelerations( const DOF & dof ) = 0;

//   //! positions
//   double & 	resultants( const DOF & dof ) = 0;

//   //! positions
//   double & 	stiffness( const DOF & dofA, const DOF & dofB ) = 0;

protected:

  Model* _model;

};

}; // namespace voom
#endif // __FESolver_h__
