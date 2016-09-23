//--------------------------------------------------------------------
//
//                  	       William Klug 
//                University of California Los Angeles
//                  (C) 2004 - All Rights Reserved
//
//--------------------------------------------------------------------
//
// Rosenbrock.cpp: Methods implementing the rosenbrock function
//
//--------------------------------------------------------------------

#include "Rosenbrock.h"
//#include <blitz/array.h>
#include <iostream>

using namespace blitz;

namespace voom
{
  // ------------
  // Constructors
  // ------------
  Rosenbrock::Rosenbrock( const int dim, const bool debug ) 
  {
    _debug = debug;
    if(_debug) std::cout << "Rosenbrock::Rosenbrock(): {" << std::endl;

    _nDOF = dim;

    _x.resize(dim);
    _grad.resize(dim);
    _x = -1.0;

    // constants
    _alpha = 105.0;

    // initialize everything else
    compute(true,true,false);

    if(_debug) printState();
  
    if(_debug) cout << "}" << std::endl;
  }

  // -------------
  // Other methods
  // -------------

  // -----------------------------------------------
  // Output to screen the current state of the model
  // -----------------------------------------------
  void Rosenbrock::printState()
  {
    Range all = Range::all();
    std::cout << "_x = " << _x << std::endl
	      << "_f = " << _f << std::endl
	      << "_grad = " << _grad << std::endl;
    
  }

  // ----------------
  // Private methods:
  // ----------------


  // ------------------------------------------------------------------
  // Compute energy density at each point on the Finite Difference grid
  // ------------------------------------------------------------------
  void Rosenbrock::compute(bool f0, bool f1, bool f2)
  {

    const int nrows = _x.rows();

    if(f0) {
      blitz::Array<double,1> temp( _x.shape() );
      blitz::Range I(0,nrows-2);
      temp(I) = 
	sqr( 1.0 - _x(I) ) + 
	_alpha*sqr( _x(I+1) - sqr(_x(I)) );
      temp(nrows-1) = 0;
      //   cout<<"Rosenbrock::_computeEnergy(){"<<endl
      //       <<"_energyDensity = "<<_energyDensity<<endl
      //       <<"}"<<endl;
      _f = sum( temp );
    }
  
    if(f1) {
  
      blitz::Array<double,1> x  = _x;
      blitz::Array<double,1> DE = _grad;
  
      DE(0) = 
	+ 2.0*( x(0) - 1.0 ) 
	- 4.0*_alpha*( x(1)-x(0)*x(0) )*x(0);

      DE(nrows-1)
	= 2.0*_alpha*( x(nrows-1)-x(nrows-2)*x(nrows-2) );

      if( nrows-2 >= 1 ) {
	Range J(1,nrows-2);
	DE(J)= 
	  + 2.0*( x(J) - 1.0 ) 
	  + 2.0*_alpha*( x(J)-x(J-1)*x(J-1) ) 
	  - 4.0*_alpha*( x(J+1)-x(J)*x(J) )*x(J);
      }
    }
    return; 
  }

  // Copy field values from nodes into arrays
 
  //! Copy positions
  void Rosenbrock::getPositions( blitz::Array< double, 1 > & x ) const {
    x = _x;
  }

  //! Copy velocities
  void Rosenbrock::getVelocities( blitz::Array< double, 1 > & v )  {
    v = 0.0;
  }
  //! Copy accelerations
  void Rosenbrock::getAccelerations( blitz::Array< double, 1 > & a )  {
    a = 0.0;
  }
  //! Copy residual forces
  void Rosenbrock::getResidual( blitz::Array< double, 1 > & r )  {
    r = -_grad;
  }
  //! Copy tangent stiffness
  void Rosenbrock::getStiffness( blitz::Array< double, 1 > & k, 
						blitz::Array< int, 2 > & ndx )  {
    k = 0.0;
    ndx = 0;
  }

  // Copy field values from arrays into nodes

  //! Copy positions
  void Rosenbrock::setPositions( const blitz::Array< double, 1 > & x ) {
    _x = x;
  }
  //! Copy velocities
  void Rosenbrock::setVelocities( const blitz::Array< double, 1 > & v ) {
  }			      
  //! Copy accelerations
  void Rosenbrock::setAccelerations( const blitz::Array< double, 1 > & a ) {
  }

}; // namespace voom  
