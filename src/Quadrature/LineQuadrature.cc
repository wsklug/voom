// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file QuadQuadrature.cc

  \brief Quadrature rule for a quadrilateral element.

*/

#include <iostream>
#include "LineQuadrature.h"

namespace voom {
  void LineQuadrature::_initialize(unsigned int order)
  {
      
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 0.0;				
      _points[0].weight = 2.0;
	
    } else if ( order == 3 ) {  // 2 points
      _points.resize(2);

      _points[0].coords = -0.577350269189626;				
      _points[1].coords =  0.577350269189626;
      _points[0].weight = _points[1].weight = 1.0;
	
    } else if ( order == 5 ) { // 3 points
      _points.resize(3);
      
      _points[0].coords = -0.774596669241483;
      _points[1].coords =  0.0;
      _points[2].coords =  0.774596669241483;
 
      _points[0].weight = 
	_points[2].weight = 5.0/9.0;

      _points[1].weight = 8.0/9.0;

    } else {
      
      std::cout << "LineQuadrature::_initialize(): No quadrature rule is implmented for order " 
		<< order << ".  Quadrature uninitialized." << std::endl;
    }
    return;
  }

  bool LineQuadrature::check(unsigned int d) const {

    // create a random polynomial of degree=d, i.e., 
    // \sum_{i+j<=d} a_{ij} s_1^i s_2^j
    srand(time(0));
    std::vector<double> a;
    for(int i=0; i<=d; i++) {
      a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.5) );
    }
    
    // Integrate the polynomial exactly using the formula from Cook's text:
    //
    //   \int_{-1}^1 \xi^i  d\xi = }
    // 
    // or with k=0
    //
    //   \int_A s_1^i s_2^j dA = 2A\frac{i!j!}{(2+i+j)!}
    //
    // For the standard quadrilateral A=1.0
    
    double I_exact=0.0;
    for(int i=0; i<=d; i++) {
      I_exact += (a[i]/(i+1))*( 1 + std::pow(-1.0,i) );
    }    

    // Integrate the polynomial P(\xi) numerically
    double I_numerical=0.0;
    for(LineQuadrature::ConstPointIterator p=this->begin(); p!=this->end(); p++) {
      double xi = p->coords(0);
      for(int i=0; i<=d; i++) {
	I_numerical += a[i]*pow(xi,i)*(p->weight);
      }
    }
    
    std::cout << "LineQuadrature::check("<<d<<"):"<<std::endl
	      << "I_exact     = " << I_exact << std::endl
	      << "I_numerical = " << I_exact << std::endl
	      << "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
	      << std::endl;

    double tol = 1.0e-8;
    if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
      std::cout << "LineQuadrature::check("<<d<<") PASSED!"
		<<std::endl;
      return true;
    }
    std::cout << "LineQuadrature::check("<<d<<") FAILED!"
	      <<std::endl;
    
    return false;
  }
} // namespace voom
