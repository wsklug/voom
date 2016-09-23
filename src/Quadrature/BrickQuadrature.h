// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
//
//----------------------------------------------------------------------

/*! 
  \file BrickQuadrature.h

  \brief Quadrature rule for a brick element.

*/

#if !defined(__BrickQuadrature_h__)
#define __BrickQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  //! Concrete Class representing 3-D Gaussian quadrature over a brick domain.
  /*! This class derives from <tt> Quadrature </tt> defining the
    method <tt>_initialize()</tt> for several quadrature orders.  In
    addition, this class provides a method <tt>check()</tt>
    which integrates exactly a random polynomial of the appropriate
    order and compares the result to that of quadrature.
  */
  class BrickQuadrature
    :public Quadrature<3>
  {

  public:

    //! default constructor
    BrickQuadrature() { _initialize(1); }
		
    //! Construct from a user-supplied quadrature order (1, 2, or 3)
    BrickQuadrature(unsigned int order) {_initialize(order);}
    
    //! assignment operator
    BrickQuadrature & operator = (const BrickQuadrature & q) {
      if( this != &q ) {
	Quadrature<3>::operator=(q);
      }
      return *this;
    }

    //! destructor
    virtual ~BrickQuadrature() {}
    bool check(unsigned int order1, unsigned int order2) const;

  private:
    //! initialize gauss points for different rules
    void _initialize(unsigned int order);  

    //! recursive funtion to compute factorial
    static int _factorial(int n) {
      int temp;
      if(n <= 1) return 1;
      temp = n * _factorial(n - 1);
      return temp;
    }
  };
	
} // namespace voom

#endif // __BrickQuadrature_h__
