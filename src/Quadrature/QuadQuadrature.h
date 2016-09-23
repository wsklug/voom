// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/04/11 04:58:17  klug
// Now log info is added directly to source file.
//
//
//----------------------------------------------------------------------

/*! 
  \file QuadQuadrature.h

  \brief Quadrature rule for a quadrilateral element.

*/

#if !defined(__QuadQuadrature_h__)
#define __QuadQuadrature_h__

#include "Quadrature.h"

namespace voom
{
  //! Concrete Class representing 2-D Gaussian quadrature over a quadrilateral domain.
  /*! This class derives from <tt> Quadrature </tt> defining the
    method <tt>_initialize()</tt> for several quadrature orders.  In
    addition, this class provides a method <tt>check()</tt>
    which integrates exactly a random polynomial of the appropriate
    order and compares the result to that of quadrature.
  */
  class QuadQuadrature
    :public Quadrature<2>
  {

  public:

    //! default constructor
    QuadQuadrature() { _initialize(1); }
		
    QuadQuadrature(unsigned int order) {_initialize(order);}
    
    //! assignment operator
    QuadQuadrature & operator = (const QuadQuadrature & q) {
      if( this != &q ) {
	Quadrature<2>::operator=(q);
      }
      return *this;
    }

    //! destructor
    virtual ~QuadQuadrature() {}
    bool check(unsigned int order) const;

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

#endif // __QuadQuadrature_h__
