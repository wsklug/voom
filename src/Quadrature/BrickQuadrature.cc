#include <iostream>
#include "BrickQuadrature.h"

namespace voom {
  void BrickQuadrature::_initialize(unsigned int order)
  {
    
    if( order == 1 ) {  // 1 point 
      _points.resize(1);

      _points[0].coords = 1.0/3.0, 1.0/3.0, 0.5;				
      _points[0].weight = 0.5;
	
    } else 

    if( order == 2) { // 2 points in the thickness and 1 point in the plane
      _points.resize(2);

      double alpha = 1.0/3.0;

      _points[0].coords = alpha, alpha, 0.211324865405187;
      _points[1].coords = alpha, alpha, 0.788675134594813;
      _points[0].weight = 0.25;
      _points[1].weight = 0.25;
 

    }  else 

    if( order == 3) { // 2 points in the thickness and 3 points in the plane
      _points.resize(6);
      double w = 1.0/12.0;

      _points[0].coords = 0.5, 0.0, 0.211324865405187;
      _points[1].coords = 0.0, 0.5, 0.211324865405187;
      _points[2].coords = 0.5, 0.5, 0.211324865405187;
      _points[3].coords = 0.5, 0.0, 0.788675134594813;
      _points[4].coords = 0.0, 0.5, 0.788675134594813;
      _points[5].coords = 0.5, 0.5, 0.788675134594813;
      _points[0].weight = w;
      _points[1].weight = w;
      _points[2].weight = w;
      _points[3].weight = w;
      _points[4].weight = w;
      _points[5].weight = w;

    } 
    return;
  }

  bool BrickQuadrature::check(unsigned int d1, unsigned int d2) const {

    BrickQuadrature quad(*this);

    // Create a random polynomial of degree=d1 in the thickness and d2 in the plane
    // In the plane //
    srand(time(0));
    std::vector<double> a;
    for(int i=0; i<=d2; i++) {
      for(int j=0; i+j<=d2; j++) {
	// a.push_back(1.0);
	a.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX - 0.33) );
      }
    }

    double I_exact_a = 0.0;
    for(int i=0, k=0; i<=d2; i++) {
      for(int j=0; i+j<=d2; j++, k++) {
	I_exact_a += a[k]*_factorial(i)*_factorial(j)/_factorial(2+i+j);
      }
    }    

    // In the thickness //
    std::vector<double> b;
    for(int i=0; i<=d1; i++) 
    {
      // b.push_back(1.0);
      b.push_back( 100.0*(static_cast<double>(rand())/RAND_MAX) - 0.5);
    }
    
    double I_exact_b = 0.0;
    for(int i=0; i<=d1; i++) 
    {
      I_exact_b += b[i]/(double(i)+1.0);
    }

    // Complete integral
    double I_exact = I_exact_a*I_exact_b;


    // Integrate the polynomial P(s) numerically
    //  \sum_p P(s_p) w_p A
    double I_numerical=0.0;
    for(BrickQuadrature::ConstPointIterator p=this->begin(); p!=this->end(); p++) {
      double s1 = p->coords(0);
      double s2 = p->coords(1);
      double s3 = p->coords(2);

      double I_numerical_a = 0.0, I_numerical_b = 0.0;

      for(int i=0, k=0; i<=d2; i++) {
	for(int j=0; i+j<=d2; j++, k++) {
	  I_numerical_a += a[k]*pow(s1,i)*pow(s2,j);
	}
      }
      
      for(int i=0; i<=d1; i++) {
	I_numerical_b += b[i]*pow(s3,i);
      }

      I_numerical += I_numerical_a*I_numerical_b*p->weight;
    }
    
    std::cout << "BrickQuadrature::check("<< d1 << " , " << d2 <<"):"<<std::endl
		<< "I_exact     = " << I_exact << std::endl
		<< "I_numerical = " << I_numerical << std::endl
		<< "Error       = " << std::abs(I_exact-I_numerical)/std::abs(I_exact) 
		<< std::endl;
	       
      double tol = 1.0e-7;
      if ( std::abs(I_numerical-I_exact) <= tol*std::abs(I_exact) ) {  
	std::cout << "BrickQuadrature::check("<< d1 << " , " << d2 <<") PASSED!"
		  <<std::endl;
	return true;
      }
      std::cout << "BrickQuadrature::check("<< d1 << " , " << d2 <<") FAILED! - FAILED! - FAILED!"
		<<std::endl;
    
      return false;
    }
  //  }
  
} // namespace voom


