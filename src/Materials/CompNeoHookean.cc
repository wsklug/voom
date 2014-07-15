// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Melissa M. Gibbons
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "VoomMath.h"
#include "CompNeoHookean.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  CompNeoHookean::CompNeoHookean(const CompNeoHookean &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _rho = Input._rho;
    _E = Input._E;
    _nu = Input._nu;
  }

	
  void CompNeoHookean::_init(double rho, double E, double nu)
  {
    _rho = rho;
    _E = E;
    _nu = nu;

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

  double CompNeoHookean::massDensity()
  {
    return _rho;
  }

  double CompNeoHookean::longitudinalWaveSpeed()
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
        
  void CompNeoHookean::updateState(bool fl0, bool fl1, bool fl2)
  {
    bool debug=false;
    
    // lame: lame constant
    // shear: shear modulus
    const double lame = _nu*_E/((1.0+_nu)*(1.0-2.0*_nu)); 
    const double shear = 0.5*_E/(1.0+_nu);
    
    // invF:  _F^-1
    Tensor3D invF;
    double jac = determinant(_F);
    
    invF(0,0) = _F(1,1)*_F(2,2) - _F(1,2)*_F(2,1);
    invF(0,1) = _F(0,2)*_F(2,1) - _F(0,1)*_F(2,2);
    invF(0,2) = _F(0,1)*_F(1,2) - _F(0,2)*_F(1,1);
    invF(1,0) = _F(1,2)*_F(2,0) - _F(1,0)*_F(2,2);
    invF(1,1) = _F(0,0)*_F(2,2) - _F(0,2)*_F(2,0);
    invF(1,2) = _F(0,2)*_F(1,0) - _F(0,0)*_F(1,2);
    invF(2,0) = _F(1,0)*_F(2,1) - _F(1,1)*_F(2,0); 
    invF(2,1) = _F(0,1)*_F(2,0) - _F(0,0)*_F(2,1);
    invF(2,2) = _F(0,0)*_F(1,1) - _F(0,1)*_F(1,0);

    invF /= jac;

    // 3-D Elasticity
    // first PK stress = (lame*ln(J) - shear)*inverse(_F) + shear*_F
    if (fl1 || fl0) {
      _P = (lame*log(jac) - shear)*(tvmet::trans(invF)) + shear*_F;
    }

    if (fl0) {
      // C:  right Cauchy-green strain tensor
      Tensor3D C;
      C = (tvmet::trans(_F))*_F;

      _W = 0.0;
      _W += 0.5*lame*(log(jac))*(log(jac)) - shear*log(jac) + 0.5*shear*(tvmet::trace(C)-3.0);
    }
  
    if (fl2) {
    }

    return;
  }


  void CompNeoHookean::ConsistencyTest()
  {
    std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
    updateState(false, true, false);
    Tensor3D PAna, PNum;  PAna = _P;                  // current 1st Piola stress
    const double eps = 1.0e-8 * max(_F);   // perturbation
    for(int i = 0; i < PAna.rows(); i ++){
      for( int j = 0; j < PAna.cols(); j++){
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
    std::cout << "Numberical value of 1st Piola stress:" << std::endl;
    std::cout << PNum << std::endl;
    std::cout << "Difference between 1st Piola stresses:" << std::endl;
    PNum -= PAna; PNum /= max(PAna);
    std::cout << PNum << std::endl;
		
    std::cout << "\n\n\n\n" << std::endl;
    return;
  }

} // namespace voom
