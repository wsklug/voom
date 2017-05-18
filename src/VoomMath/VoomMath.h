#ifndef _VOOM_MATH_
#define _VOOM_MATH_
 
#include "voom.h"
#include <Eigen/Geometry>

namespace voom
{

  inline double sqr(double a) {return (a*a);}

//   inline double norm2(const blitz::TinyVector<double,3> & v) { 
//     return std::sqrt( v(0)*v(0) + v(1)*v(1) + v(2)*v(2) );
//   }

//   inline double dot(const blitz::TinyVector<double,3> & u, const blitz::TinyVector<double,3> & v) {
//     return u(0)*v(0) + u(1)*v(1) + u(2)*v(2);
//   }
  
  double norm2(const Array1D  & v);

  void invert(const Tensor3D & a, Tensor3D & b);

  void invert(const Tensor2D & a, Tensor2D & b);

  double determinant(const Tensor2D & F);

  double determinant(const Tensor3D & F);

  void tensorProduct(const Vector3D u, const Vector3D v, Tensor3D & T);

  void tensorProduct(const Vector2D u, const Vector2D v, Tensor2D & T);

  Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out);

  void TestFind3DAffineTransform();

};

#endif // _VOOM_MATH_
