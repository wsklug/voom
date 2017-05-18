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
  \file QuadQuadrature.cc

  \brief Quadrature rule for a quadrilateral element.

*/

#include <iostream>
#include "QuadQuadrature.h"

namespace voom {
  void QuadQuadrature::_initialize(unsigned int order)
  {
      
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 0.0, 0.0;				
      _points[0].weight = 4.0;
	
    } else if ( order == 2 ) {  // 4 points
      _points.resize(4);

      _points[0].coords = -0.577350269, -0.577350269;				
      _points[1].coords = 0.577350269, -0.577350269;				
      _points[2].coords = 0.577350269, 0.577350269;
      _points[3].coords = -0.577350269, 0.577350269;
      _points[0].weight = 
	_points[1].weight =
	_points[2].weight =
	_points[3].weight = 1.0;
	
    } else if ( order == 3 ) { // 9 points
      _points.resize(9);
      
      
      _points[0].coords = -0.244948974,-0.244948974;
      _points[1].coords =  0.0        ,-0.244948974;
      _points[2].coords =  0.244948974,-0.244948974;
      _points[3].coords = -0.244948974, 0.0;
      _points[4].coords =  0.0        , 0.0;
      _points[5].coords =  0.244948974, 0.0;
      _points[6].coords = -0.244948974, 0.244948974;
      _points[7].coords =  0.0        , 0.244948974;
      _points[8].coords =  0.244948974, 0.244948974;

      _points[4].weight = 64.0/81.0;
 
      _points[0].weight =
	_points[2].weight = 
	_points[6].weight = 
	_points[8].weight = 25.0/81.0;

      _points[1].weight =
	_points[3].weight = 
	_points[5].weight = 
	_points[7].weight = 40.0/81.0;

    } else {
      
      std::cout << "QuadQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    }
    return;
  }

  bool QuadQuadrature::check(unsigned int d) const {

    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(0));
    std::vector<double> a;
    for(int i=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++) {
	a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.33) );
      }
    }
    
    // Integrate the polynomial exactly using the formula from Cook's text:
    //
    //   \int_A s_1^i s_2^j (1-s_1-s_2)^k dA = 2A\frac{i!j!m!}{(2+i+j+k)!}
    // 
    // or with k=0
    //
    //   \int_A s_1^i s_2^j dA = 2A\frac{i!j!}{(2+i+j)!}
    //
    // For the standard quadrilateral A=1.0
    
    double I_exact=0.0;
    for(int i=0, k=0; i<=d; i++) {
      for(int j=0; i+j<=d; j++, k++) {
	I_exact += a[k]*_factorial(i)*_factorial(j)/_factorial(2+i+j);
      }
    }    

    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double I_numerical=0.0;
    for(QuadQuadrature::ConstPointIterator p=this->begin(); p!=this->end(); p++) {
      double s1 = p->coords(0);
      double s2 = p->coords(1);
      for(int i=0, k=0; i<=d; i++) {
	for(int j=0; i+j<=d; j++, k++) {
	  I_numerical += a[k]*pow(s1,i)*pow(s2,j)*(p->weight);
	}
      }
    }
    I_numerical *= 0.5;
    
    std::cout << "QuadQuadrature::check("<<d<<"):"<<std::endl
	      << "I_exact     = " << I_exact << std::endl
	      << "I_numerical = " << I_exact << std::endl
	      << "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
	      << std::endl;

    double tol = 1.0e-8;
    if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
      std::cout << "QuadQuadrature::check("<<d<<") PASSED!"
		<<std::endl;
      return true;
    }
    std::cout << "QuadQuadrature::check("<<d<<") FAILED!"
	      <<std::endl;
    
    return false;
  }
} // namespace voom
