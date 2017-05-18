#include "VoomMath.h"

namespace voom {

  double norm2(Array1D & v)  {
    double result = 0.0;
    for( int i = 0; i < v.size(); i++){
      result += sqr( v(i) );
    }
    return std::sqrt(result);
  }
      	
  double determinant(const Tensor2D & F)  {
    return ( F(0,0)*F(1,1) - F(0,1)*F(1,0) );
  }

  double determinant(const Tensor3D & F)  {
    return (double)( F(0,0)*F(1,1)*F(2,2) +
		     F(1,0)*F(2,1)*F(0,2) +
		     F(0,1)*F(1,2)*F(2,0) -
		     F(0,2)*F(1,1)*F(2,0) -
		     F(0,1)*F(1,0)*F(2,2) -
		     F(1,2)*F(2,1)*F(0,0) );
  }

	
  void invert(const Tensor2D & a, Tensor2D & b)  {
    const double det = determinant(a);

    b(0,0) =  a(1,1)/det;
    b(0,1) = -a(0,1)/det;
    b(1,0) = -a(1,0)/det;
    b(1,1) =  a(0,0)/det;
  
    return;
  }

  void invert(const Tensor3D & a, Tensor3D & b)  {
    const double det = determinant(a);

    b(0,0) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2) )/det;
    b(1,0) = ( a(2,0)*a(1,2) - a(1,0)*a(2,2) )/det;
    b(2,0) = ( a(1,0)*a(2,1) - a(1,1)*a(2,0) )/det;
    b(0,1) = ( a(2,1)*a(0,2) - a(0,1)*a(2,2) )/det;
    b(1,1) = ( a(0,0)*a(2,2) - a(2,0)*a(0,2) )/det;
    b(2,1) = ( a(0,1)*a(2,0) - a(0,0)*a(2,1) )/det;
    b(0,2) = ( a(0,1)*a(1,2) - a(1,1)*a(0,2) )/det;
    b(1,2) = ( a(0,2)*a(1,0) - a(0,0)*a(1,2) )/det;
    b(2,2) = ( a(0,0)*a(1,1) - a(1,0)*a(0,1) )/det;
  
    return;
  }

  void tensorProduct(const Vector3D u, const Vector3D v, Tensor3D & T) {    
    for(int i=0; i<3; i++) 
      for(int j=0; j<3; j++) 
	T(i,j) = u(i)*v(j);
    return;
  } 

  void tensorProduct(const Vector2D u, const Vector2D v, Tensor2D & T) {    
    for(int i=0; i<2; i++) 
      for(int j=0; j<2; j++) 
	T(i,j) = u(i)*v(j);
    return;
  } 

	
};
