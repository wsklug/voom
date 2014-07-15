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
// Revision 1.3  2005/05/23 17:39:11  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

#include <iostream>
#include "VoomMath.h"
#include "StVenant.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  StVenant::StVenant(const StVenant &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _rho = Input._rho;
    _E = Input._E;
    _nu = Input._nu;
  }

	
  void StVenant::_init(double rho, double E, double nu)
  {
    _rho = rho;
    _E = E;
    _nu = nu;

    //std::cout << "_nu = " << _nu << std::endl;
    //
    // deformation gradient tensor
    _F =
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;
    //
    // strain energy density
    _W = 0.0;
    //
    // first Piola stress tensor
    _P = 
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0,
      0.0, 0.0, 0.0;
  }

  // Accessors/mutators

  double StVenant::massDensity()
  {
    return _rho;
  }

  double StVenant::longitudinalWaveSpeed()
  {
    // lame: lame constant
    // shear: shear modulus
    const double lame = _nu*_E/((1.0+_nu)*(1.0-2.0*_nu)); 
    const double shear = 0.5*_E/(1.0+_nu);
    //                       ______________________
    //                      /
    //                     /  lame + 2.0 * shear
    //  P wave speed  ==  /  ____________________
    //                   /
    //                 \/            rho
    //
    return sqrt((lame + 2.0 * shear)/_rho);
  }

  // Operators


  // General methods
        
  void StVenant::updateState(bool fl0, bool fl1, bool fl2)
  {
    bool debug=false;
    //
    // E: Green Strain Tensor
    // S:       Stress Tensor
    Tensor3D E, S;
    E =	1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0;

    //E+= (tvmet::trans(_F))*_F;
    E = (tvmet::trans(_F))*_F - E;		
    E /= 2.0;  
  
    if(debug){
      std::cout << "E(I,J) = " << std::endl;
      std::cout << E << std::endl;
    }

    //
    // lame: lame constant
    // shear: shear modulus
    const double lame = _nu*_E/((1.0+_nu)*(1.0-2.0*_nu)); 
    const double shear = 0.5*_E/(1.0+_nu);
    const double trace = tvmet::trace(E);

    //
    // 3-D Elasticity
    // stress = lame * trace(strain) * Identity + 2*shear*strain
    if (fl1 || fl0) {
      S = 2.0*shear*E;
      for(int J=0; J<3; J++) {
	S(J,J) += lame*trace;
      }
		
      if(debug){
	std::cout << "S(J,I) = " << std::endl;
	for(int J=0; J<3; J++) {
	  std::cout << "\t";
	  for(int I=0; I<3; I++) {
	    std::cout << S(J,I) << "  ";
	  }
	  std::cout << std::endl;
	}
	std::cout << std::endl;
      }
    }

    if (fl0) {
      _W = 0.0;
      //_W += 0.5*tvmet::trace(S*E);
      _W += 0.5*tvmet::sum(tvmet::diag(S*E));
    }

    if (fl1) {
      _P = _F*S;

      double jac = determinant(_F);
      _C = 1.0/jac*_P*(tvmet::trans(_F));

      _vMises = 0.7071068*sqrt( ( _C(0,0) - _C(1,1) )*( _C(0,0) - _C(1,1) ) + 
				( _C(1,1) - _C(2,2) )*( _C(1,1) - _C(2,2) ) + 
				( _C(2,2) - _C(0,0) )*( _C(2,2) - _C(0,0) ) + 
				6.0*( _C(0,1)*_C(0,1) + _C(0,2)*_C(0,2) + _C(1,2)*_C(1,2) ) );
    }
  
    if (fl2) {
      //     for( int L=0; L<3; L++)
      //       for( int k=0; k<3; k++)
      // 	for( int J=0; J<3; J++)
      // 	  for( int i=0; i<3; i++)
      // 	    DDw(L,k,J,i) 
      // 	      = lame*F(3*J+i)*F(3*L+k) 
      // 	      + shear*F(3*L+i)*F(3*J+k);
    
      //     for( int k=0; k<3; k++) 
      //       for( int i=0; i<3; i++) {
      // 	double FiK_FkK = 0.0;
      // 	for( int K=0; K<3; K++) FiK_FkK += F(3*K+i)*F(3*K+k);
      // 	for( int J=0; J<3; J++) DDw(J,k,J,i) += shear*FiK_FkK;
      //       }
    
      //     for( int J=0; J<3; J++)
      //       for( int i=0; i<3; i++) 
      // 	DDw(J,i,J,i) += lame*trace;
    
      //     for( int L=0; L<3; L++)
      //       for( int J=0; J<3; J++)
      // 	for( int i=0; i<3; i++) 
      // 	  DDw(L,i,J,i) += 2.0*shear*E(L,J);
    }
    return;
  }


  void StVenant::ConsistencyTest()
  {
    std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
    updateState(false, true, false);
    Tensor3D PAna, PNum;  PAna = _P;                  // current 1st Piola stress
    const double eps = 1.0e-8 * max(_F);   // perturbation
    for(int i = 0; i < 3; i ++){
      for( int j = 0; j < 3; j++){
	_F(i,j) += eps;
	updateState(true, false, false);
	double W = _W;
	_F(i,j) -= 2*eps;
	updateState(true, false, false);
	W -= _W;
	_F(i,j) += eps;       // restore value
	PNum(i,j) = W/2.0/eps;
      }
    }
    std::cout << "Analytical value of 1st Piola stress:" << std::endl;
    std::cout << PAna << std::endl;
    std::cout << "Numerical value of 1st Piola stress:" << std::endl;
    std::cout << PNum << std::endl;
    std::cout << "Difference between 1st Piola stresses:" << std::endl;
    PNum -= PAna; PNum /= max(PAna);
    std::cout << PNum << std::endl;
		
    std::cout << "\n\n\n\n" << std::endl;
    return;
  }

  // 	double StVenant::_determinant(const Tensor3D & F)
  // 	{
  // 		double J;

  // 		J = F(0,0)*F(1,1)*F(2,2) + F(1,0)*F(2,1)*F(0,2) + F(0,1)*F(1,2)*F(2,0) 
  // 			- F(0,2)*F(1,1)*F(2,0) - F(0,1)*F(1,0)*F(2,2) - F(1,2)*F(2,1)*F(0,0);

  // 		return J;
  // 	}

  // 	void StVenant::_invert(const Tensor3D & a, Tensor3D & b)
  // 	{
  // 		double det = _determinant(a);

  // 		b(0,0) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2) )/det;
  // 		b(1,0) = ( a(2,0)*a(1,2) - a(1,0)*a(2,2) )/det;
  // 		b(2,0) = ( a(1,0)*a(2,1) - a(1,1)*a(2,0) )/det;
  // 		b(0,1) = ( a(2,1)*a(0,2) - a(0,1)*a(2,2) )/det;
  // 		b(1,1) = ( a(0,0)*a(2,2) - a(2,0)*a(0,2) )/det;
  // 		b(2,1) = ( a(0,1)*a(2,0) - a(0,0)*a(2,1) )/det;
  // 		b(0,2) = ( a(0,1)*a(1,2) - a(1,1)*a(0,2) )/det;
  // 		b(1,2) = ( a(0,2)*a(1,0) - a(0,0)*a(1,2) )/det;
  // 		b(2,2) = ( a(0,0)*a(1,1) - a(1,0)*a(0,1) )/det;
  
  // 		return;
  // 	}

} // namespace voom
