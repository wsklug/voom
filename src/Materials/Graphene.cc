// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2012 All Rights Reserved
//
//----------------------------------------------------------------------

#include <iostream>
#include "Graphene.h"

namespace voom 
{
  void Graphene::updateState(bool f0, bool f1, bool f2)
  {
    // Elastic coefficients from Wei, et al., PRB 80: 205407 (2009)
    //
    // units = pN/nm = 1e-3 N/m
    //
    // SOEC
    //
    // linear moduli are temporarily reduced by 20% to see effect
    //
    const double C11=300.0e3; //358.1e3;
    const double C12=50.0e3;  //60.4e3;
 
    // TOEC
    const double C111=-2817.0e3;
    const double C112=-337.1e3;
    const double C222=-2693.3e3;

    // FOEC
    const double C1111=13416.2e3;
    const double C1112=759.0e3;
    const double C1122=2582.8e3;
    const double C2222=10358.9e3;

    // FFOEC
    const double C11111=-31383.8e3;
    const double C11112=-88.4e3;
    const double C11122=-12960.5e3;
    const double C12222=-13046.6e3;
    const double C22222=-33446.7e3;

    // // TOEC
    // const double C111=0.0;
    // const double C112=0.0;
    // const double C222=0.0;

    // // FOEC
    // const double C1111=0.0;
    // const double C1112=0.0;
    // const double C1122=0.0;
    // const double C2222=0.0;

    // // FFOEC
    // const double C11111=0.0;
    // const double C11112=0.0;
    // const double C11122=0.0;
    // const double C12222=0.0;
    // const double C22222=0.0;

    // basis vectors

    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;

    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();

    const Tensor2D& metric = _deformedGeometry.metricTensor();
    const Tensor2D& refMetric = _referenceGeometry.metricTensor();
    const Tensor2D& refMetricInv = _referenceGeometry.metricTensorInverse();

    _strain = 0.5*(metric-refMetric);

    // Tangent Cartesian basis
    BasisVectors e;

 //   e(0) = refBasis(0)/norm2(refBasis(0));
 //   e(1) = refDual(1)/norm2(refDual(1));
    // For now, use x and y axes in the "Lab frame"
    // Should change this at some point to allow for a more general definition of anisotropy.
    // Maybe pass in the special directions?

    // e(0) = 1.0, 0.0, 0.0;
    // e(1) = 0.0, 1.0, 0.0;
    double c = cos(_theta);
    double s = sin(_theta);

    e(0) =  c, s, 0.0;
    e(1) = -s, c, 0.0;

    double eta1 = 0.0;
    double eta2 = 0.0;
    double eta3 = 0.0;

    // Transformation matrix between Cartesian and curvilinear bases
    Tensor2D T(0.0);
    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	T(alpha,beta) = dot( refDual(alpha), e(beta) );
      }
    }
    for(int alpha=0; alpha<2; alpha++) {
      for(int beta=0; beta<2; beta++) {
	eta1 += _strain(alpha,beta)*T(alpha,0)*T(beta,0);
	eta2 += _strain(alpha,beta)*T(alpha,1)*T(beta,1);
	eta3 += _strain(alpha,beta)*(T(alpha,0)*T(beta,1)+T(alpha,1)*T(beta,0));
      }
    }
    
    // Tensor2D strainDual(0.0);
    // strainDual = refMetricInv*_strain*tvmet::trans(refMetricInv);

    // double traceStrain = 0.0;
    // double strainSquared = 0.0;
    // for(int alpha=0; alpha<2; alpha++) {
    //   for(int beta=0; beta<2; beta++) {
    // 	traceStrain   += refMetricInv(alpha,beta)*_strain(alpha,beta);
    // 	strainSquared += strainDual(alpha,beta)*_strain(alpha,beta);
    //   }
    // }

    // Vector3D xhat; xhat = 1.0, 0.0, 0.0;
    // Vector3D yhat; yhat = 0.0, 1.0, 0.0;

    // Vector2D ex(0.0);
    // for(int alpha=0; alpha<2; alpha++) ex(alpha) = dot(xhat, refDual(alpha));

    // Vector2D ey(0.0);
    // for(int alpha=0; alpha<2; alpha++) ey(alpha) = dot(yhat, refDual(alpha));

    // double epsxx = 0.0;
    // double epsyy = 0.0;
    // double epsxy = 0.0;

    // for(int alpha=0; alpha<2; alpha++) {
    //   for(int beta=0; beta<2; beta++) {
    // 	epsxx += _strain(alpha,beta)*ex(alpha)*ex(beta);
    // 	epsyy += _strain(alpha,beta)*ey(alpha)*ey(beta);
    // 	epsxy += _strain(alpha,beta)*ex(alpha)*ey(beta)
    // 	      +  _strain(alpha,beta)*ey(alpha)*ex(beta);
    //   }
    // }
    // epsxy *= 0.5;

    // // contravariant components of 2nd P-K stress
    // _stress = _lambda*traceStrain*refMetricInv + 2.0*_mu*strainDual;

    // // 3rd order terms
    // double DuDepsxx 
    //   = 0.5*_C111*epsxx*epsxx 
    //   + _C112*epsxx*epsyy
    //   + 0.5*(_C111-_C222+_C112)*epsyy*epsyy
    //   + 0.5*(3.0*_C222-2.0*_C111-_C112)*epsxy*epsxy;

    // double DuDepsyy 
    //   = 0.5*_C222*epsyy*epsyy 
    //   + 0.5*_C112*epsxx*epsxx
    //   + (_C111-_C222+_C112)*epsxx*epsyy
    //   + 0.5*(2.0*_C111-_C222-_C112)*epsxy*epsxy;

    // double DuDepsxy 
    //   = (3.0*_C222-2.0*_C111-_C112)*epsxx*epsxy
    //   + (2.0*_C111-_C222-_C112)*epsyy*epsxy;

    // for(int alpha=0; alpha<2; alpha++) {
    //   for(int beta=0; beta<2; beta++) {
    // 	_stress(alpha,beta) += 
    // 	  DuDepsxx*ex(alpha)*ex(beta) +
    // 	  DuDepsyy*ey(alpha)*ey(beta) +
    // 	  0.5*DuDepsxy*( ex(alpha)*ey(beta)+ey(alpha)*ex(beta) );
    //   }
    // }

    double jacobian = _referenceGeometry.metric()/_deformedGeometry.metric();

    if( f0 ) {
      // compute strain energy
      // _W = (0.5*_lambda*traceStrain*traceStrain + _mu*strainSquared)*jacobian;

      // // 3rd-order terms
      // _W += ( _C111*pow(epsxx,3)/6.0 +  _C222*pow(epsyy,3)/6.0
      // 	+ 0.5*_C112*epsxx*epsxx*epsyy
      // 	+ 0.5*(_C111-_C222+_C112)*epsxx*epsyy*epsyy
      // 	+ 0.5*(3.0*_C222-2.0*_C111-_C112)*epsxx*epsxy*epsxy
      // 	      + 0.5*(2.0*_C111-_C222-_C112)*epsyy*epsxy*epsxy ) * jacobian;

      _W = (C11111*pow(eta1,5) + 
     5*C11112*pow(eta1,4)*eta2 + 
     10*C11122*pow(eta1,3)*pow(eta2,2) + 
     5*(C11111 + 3*C11112 + 2*C11122 - 3*C12222 - 
        C22222)*pow(eta1,2)*pow(eta2,3) + 
     5*C12222*eta1*pow(eta2,4) + 
     C22222*pow(eta2,5) - 
     ((4*C11111 + 5*C11112 - 9*C22222)*pow(eta1,3)*
        pow(eta3,2))/4. - 
     ((13*C11111 + 30*C11112 + 20*C11122 - 
          45*C12222 - 18*C22222)*pow(eta1,2)*eta2*
        pow(eta3,2))/4. + 
     ((8*C11111 + 15*C11112 - 20*C11122 - 3*C22222)*
        eta1*pow(eta2,2)*pow(eta3,2))/4. + 
     ((9*C11111 - 5*C12222 - 4*C22222)*pow(eta2,3)*
        pow(eta3,2))/4. + 
     ((11*C11111 + 30*C11112 + 10*C11122 - 
          45*C12222 - 6*C22222)*eta1*pow(eta3,4))/
      16. - ((C11111 + 30*C11112 - 10*C11122 - 
          15*C12222 - 6*C22222)*eta2*pow(eta3,4))/
      16. + 60*(C11*pow(eta1,2) + 
        2*C12*eta1*eta2 + C11*pow(eta2,2) + 
        ((C11 - C12)*pow(eta3,2))/2.) + 
     20*(C111*pow(eta1,3) + 
        3*C112*pow(eta1,2)*eta2 + 
        3*(C111 + C112 - C222)*eta1*pow(eta2,2) + 
        C222*pow(eta2,3) - 
        (3*(2*C111 + C112 - 3*C222)*eta1*
           pow(eta3,2))/4. + 
        (3*(2*C111 - C112 - C222)*eta2*
           pow(eta3,2))/4.) + 
     5*(C1111*pow(eta1,4) + 
        4*C1112*pow(eta1,3)*eta2 + 
        6*C1122*pow(eta1,2)*pow(eta2,2) + 
        2*(C1111 + 2*C1112 - C2222)*eta1*
         pow(eta2,3) + C2222*pow(eta2,4) - 
        ((5*C1111 + 4*C1112 - 9*C2222)*pow(eta1,2)*
           pow(eta3,2))/4. + 
        (C1111 + 2*C1112 - 3*C1122)*eta1*eta2*
         pow(eta3,2) + 
        ((7*C1111 - 4*C1112 - 3*C2222)*pow(eta2,2)*
           pow(eta3,2))/4. - 
        ((C1111 + 8*C1112 - 3*(2*C1122 + C2222))*
	 pow(eta3,4))/16.))/120.;

      _W *= jacobian;

    }
		
    if( f1 ) {

      double S11 = (5*C11111*pow(eta1,4) + 
     20*C11112*pow(eta1,3)*eta2 + 
     30*C11122*pow(eta1,2)*pow(eta2,2) + 
     10*(C11111 + 3*C11112 + 2*C11122 - 3*C12222 - 
        C22222)*eta1*pow(eta2,3) + 
     5*C12222*pow(eta2,4) + 
     120*(C11*eta1 + C12*eta2) - 
     (3*(4*C11111 + 5*C11112 - 9*C22222)*
        pow(eta1,2)*pow(eta3,2))/4. - 
     ((13*C11111 + 30*C11112 + 20*C11122 - 
          45*C12222 - 18*C22222)*eta1*eta2*
        pow(eta3,2))/2. + 
     ((8*C11111 + 15*C11112 - 20*C11122 - 3*C22222)*
        pow(eta2,2)*pow(eta3,2))/4. + 
     ((11*C11111 + 30*C11112 + 10*C11122 - 
          45*C12222 - 6*C22222)*pow(eta3,4))/16. + 
     20*(3*C111*pow(eta1,2) + 6*C112*eta1*eta2 + 
        3*(C111 + C112 - C222)*pow(eta2,2) - 
        (3*(2*C111 + C112 - 3*C222)*pow(eta3,2))/4.
        ) + 5*(4*C1111*pow(eta1,3) + 
        12*C1112*pow(eta1,2)*eta2 + 
        12*C1122*eta1*pow(eta2,2) + 
        2*(C1111 + 2*C1112 - C2222)*pow(eta2,3) - 
        ((5*C1111 + 4*C1112 - 9*C2222)*eta1*
           pow(eta3,2))/2. + 
        (C1111 + 2*C1112 - 3*C1122)*eta2*
	       pow(eta3,2)))/120.;

      double S22 = (5*C11112*pow(eta1,4) + 
     20*C11122*pow(eta1,3)*eta2 + 
     15*(C11111 + 3*C11112 + 2*C11122 - 3*C12222 - 
        C22222)*pow(eta1,2)*pow(eta2,2) + 
     20*C12222*eta1*pow(eta2,3) + 
     5*C22222*pow(eta2,4) + 
     120*(C12*eta1 + C11*eta2) - 
     ((13*C11111 + 30*C11112 + 20*C11122 - 
          45*C12222 - 18*C22222)*pow(eta1,2)*
        pow(eta3,2))/4. + 
     ((8*C11111 + 15*C11112 - 20*C11122 - 3*C22222)*
        eta1*eta2*pow(eta3,2))/2. + 
     (3*(9*C11111 - 5*C12222 - 4*C22222)*
        pow(eta2,2)*pow(eta3,2))/4. - 
     ((C11111 + 30*C11112 - 10*C11122 - 15*C12222 - 
          6*C22222)*pow(eta3,4))/16. + 
     20*(3*C112*pow(eta1,2) + 
        6*(C111 + C112 - C222)*eta1*eta2 + 
        3*C222*pow(eta2,2) + 
        (3*(2*C111 - C112 - C222)*pow(eta3,2))/4.)\
      + 5*(4*C1112*pow(eta1,3) + 
        12*C1122*pow(eta1,2)*eta2 + 
        6*(C1111 + 2*C1112 - C2222)*eta1*
         pow(eta2,2) + 4*C2222*pow(eta2,3) + 
        (C1111 + 2*C1112 - 3*C1122)*eta1*
         pow(eta3,2) + 
        ((7*C1111 - 4*C1112 - 3*C2222)*eta2*
	 pow(eta3,2))/2.))/120.;

      double S12 = (eta3*(240*(C11 - C12) - 
       2*(4*C11111 + 5*C11112 - 9*C22222)*
        pow(eta1,3) - 
       2*(13*C11111 + 30*C11112 + 20*C11122 - 
          45*C12222 - 18*C22222)*pow(eta1,2)*eta2\
        + 2*(8*C11111 + 15*C11112 - 20*C11122 - 
          3*C22222)*eta1*pow(eta2,2) + 
       2*(9*C11111 - 5*C12222 - 4*C22222)*
        pow(eta2,3) - 
       120*(2*C111*(eta1 - eta2) + 
          C222*(-3*eta1 + eta2) + C112*(eta1 + eta2))
         + (11*C11111 + 30*C11112 + 10*C11122 - 
          45*C12222 - 6*C22222)*eta1*pow(eta3,2) + 
       (-C11111 - 30*C11112 + 10*C11122 + 
          15*C12222 + 6*C22222)*eta2*pow(eta3,2) + 
       5*(-2*(5*C1111 + 4*C1112 - 9*C2222)*
           pow(eta1,2) + 
          8*(C1111 + 2*C1112 - 3*C1122)*eta1*eta2 + 
          2*(7*C1111 - 4*C1112 - 3*C2222)*
           pow(eta2,2) - 
          (C1111 + 8*C1112 - 3*(2*C1122 + C2222))*
	  pow(eta3,2))))/480.;

      _stress = Tensor2D(0.0);

      for(int alpha=0; alpha<2; alpha++) {
	for(int beta=0; beta<2; beta++) {
	  _stress(alpha,beta) += 
	    S11*T(alpha,0)*T(beta,0) +
	    S22*T(alpha,1)*T(beta,1) +
	    S12*(T(alpha,0)*T(beta,1)+T(alpha,1)*T(beta,0));
	}
      }
      // stress resultants
      for(int alpha=0; alpha<2; alpha++) {
	_n(alpha) = 0.0, 0.0, 0.0;
	_m(alpha) = 0.0, 0.0, 0.0;
	for(int beta=0; beta<2; beta++) {
	  _n(alpha) += _stress(alpha,beta) * basis(beta) * jacobian;
	}
      }

      _n(2) = 0.0, 0.0, 0.0;
    }

    return;
  }
}
