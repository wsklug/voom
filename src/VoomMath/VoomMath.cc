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

  // some routines borrowed from the Java Matrix library JAMA
  // needed for eigenDecomposition()


  double hypot2(double x, double y) {
    return sqrt(x*x+y*y);
  }

  // Symmetric Householder reduction to tridiagonal form.

  void tred2(double V[3][3], double d[3], double e[3]) {

    const int n=3;

    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--) {

      // Scale to avoid under/overflow.

      double scale = 0.0;
      double h = 0.0;
      for (int k = 0; k < i; k++) {
	scale = scale + fabs(d[k]);
      }
      if (scale == 0.0) {
	//if (scale < 1.0e-14) {
	e[i] = d[i-1];
	for (int j = 0; j < i; j++) {
	  d[j] = V[i-1][j];
	  V[i][j] = 0.0;
	  V[j][i] = 0.0;
	}
      } else {

	// Generate Householder vector.

	for (int k = 0; k < i; k++) {
	  d[k] /= scale;
	  h += d[k] * d[k];
	}
	double f = d[i-1];
	double g = sqrt(h);
	if (f > 0) {
	  g = -g;
	}
	e[i] = scale * g;
	h = h - f * g;
	d[i-1] = f - g;
	for (int j = 0; j < i; j++) {
	  e[j] = 0.0;
	}

	// Apply similarity transformation to remaining columns.

	for (int j = 0; j < i; j++) {
	  f = d[j];
	  V[j][i] = f;
	  g = e[j] + V[j][j] * f;
	  for (int k = j+1; k <= i-1; k++) {
	    g += V[k][j] * d[k];
	    e[k] += V[k][j] * f;
	  }
	  e[j] = g;
	}
	f = 0.0;
	for (int j = 0; j < i; j++) {
	  e[j] /= h;
	  f += e[j] * d[j];
	}
	double hh = f / (h + h);
	for (int j = 0; j < i; j++) {
	  e[j] -= hh * d[j];
	}
	for (int j = 0; j < i; j++) {
	  f = d[j];
	  g = e[j];
	  for (int k = j; k <= i-1; k++) {
	    V[k][j] -= (f * e[k] + g * d[k]);
	  }
	  d[j] = V[i-1][j];
	  V[i][j] = 0.0;
	}
      }
      d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
      V[n-1][i] = V[i][i];
      V[i][i] = 1.0;
      double h = d[i+1];
      if (h != 0.0) {
	for (int k = 0; k <= i; k++) {
	  d[k] = V[k][i+1] / h;
	}
	for (int j = 0; j <= i; j++) {
	  double g = 0.0;
	  for (int k = 0; k <= i; k++) {
	    g += V[k][i+1] * V[k][j];
	  }
	  for (int k = 0; k <= i; k++) {
	    V[k][j] -= g * d[k];
	  }
	}
      }
      for (int k = 0; k <= i; k++) {
	V[k][i+1] = 0.0;
      }
    }
    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
      V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
  } 

  // Symmetric tridiagonal QL algorithm.

  static void tql2(double V[3][3], double d[3], double e[3]) {

    const int n=3;

    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (int i = 1; i < n; i++) {
      e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (int l = 0; l < n; l++) {

      // Find small subdiagonal element

      tst1 = std::max(tst1,fabs(d[l]) + fabs(e[l]));
      int m = l;
      while (m < n) {
	if (fabs(e[m]) <= eps*tst1) {
	  break;
	}
	m++;
      }

      // If m == l, d[l] is an eigenvalue,
      // otherwise, iterate.

      if (m > l) {
	int iter = 0;
	do {
	  iter = iter + 1;  // (Could check iteration count here.)

	  // Compute implicit shift

	  double g = d[l];
	  double p = (d[l+1] - g) / (2.0 * e[l]);
	  double r = hypot2(p,1.0);
	  if (p < 0) {
	    r = -r;
	  }
	  d[l] = e[l] / (p + r);
	  d[l+1] = e[l] * (p + r);
	  double dl1 = d[l+1];
	  double h = g - d[l];
	  for (int i = l+2; i < n; i++) {
	    d[i] -= h;
	  }
	  f = f + h;

	  // Implicit QL transformation.

	  p = d[m];
	  double c = 1.0;
	  double c2 = c;
	  double c3 = c;
	  double el1 = e[l+1];
	  double s = 0.0;
	  double s2 = 0.0;
	  for (int i = m-1; i >= l; i--) {
	    c3 = c2;
	    c2 = c;
	    s2 = s;
	    g = c * e[i];
	    h = c * p;
	    r = hypot2(p,e[i]);
	    e[i+1] = s * r;
	    s = e[i] / r;
	    c = p / r;
	    p = c * d[i] - s * g;
	    d[i+1] = h + s * (c * g + s * d[i]);

	    // Accumulate transformation.

	    for (int k = 0; k < n; k++) {
	      h = V[k][i+1];
	      V[k][i+1] = s * V[k][i] + c * h;
	      V[k][i] = c * V[k][i] - s * h;
	    }
	  }
	  p = -s * s2 * c3 * el1 * e[l] / dl1;
	  e[l] = s * p;
	  d[l] = c * p;

	  // Check for convergence.

	} while (fabs(e[l]) > eps*tst1);
      }
      d[l] = d[l] + f;
      e[l] = 0.0;
    }
  
    // Sort eigenvalues and corresponding vectors.

    for (int i = 0; i < n-1; i++) {
      int k = i;
      double p = d[i];
      for (int j = i+1; j < n; j++) {
	if (d[j] < p) {
	  k = j;
	  p = d[j];
	}
      }
      if (k != i) {
	d[k] = d[i];
	d[i] = p;
	for (int j = 0; j < n; j++) {
	  p = V[j][i];
	  V[j][i] = V[j][k];
	  V[j][k] = p;
	}
      }
    }
  }


  void eigenDecomposition( const Tensor3D & A, 
			   Tensor3D & eVecs, 
			   Vector3D & eVals   ) {
    const int n=3;
    double e[n];
    double V[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	V[i][j] = A(i,j);
      }
    }

    double *d = &(eVals[0]);
    tred2(V, d, e);

    std::cout << "d=" << std::endl;
    for (int i = 0; i < n; i++) {
      std::cout << d[i] << " " ;
    }
    std::cout << std::endl;
    std::cout << "e=" << std::endl;
    for (int i = 0; i < n; i++) {
      std::cout << e[i] << " " ;
    }
    std::cout << std::endl;
    std::cout << "V=" << std::endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	std::cout << V[i][j] << " " ;
      }
    }
    std::cout << std::endl;

    tql2(V, d, e);

    std::cout << "d=" << std::endl;
    for (int i = 0; i < n; i++) {
      std::cout << d[i] << " " ;
    }
    std::cout << std::endl;
    std::cout << "V=" << std::endl;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	std::cout << V[i][j] << " " ;
      }
    }
    std::cout << std::endl;

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	eVecs(i,j) = V[i][j];
      }
    }    
  }
	
};
