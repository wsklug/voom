// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  Ankush Aggarwal, William S. Klug
//                University of California Los Angeles
//                   (C) 2011 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include <iostream>
#include "VoomMath.h"
#include "HomogMP.h"

//using namespace std;
namespace voom {

  // Construction/Destruction
  HomogMP::HomogMP(const HomogMP &Input)
  {
    if (this == &Input) return;

    _F = Input._F;

    _C10 = Input._C10;
    _C01 = Input._C01;
    _D1 = Input._D1;

    _I1MP = Input._I1MP;
    _I2MP = Input._I2MP;
    _JMP = Input._JMP;
  }

	
  void HomogMP::_init(double C10, double C01, double D1)
  {
    _C10 = C10;
    _C01 = C01;
    _D1 = D1;

    _I1MP = 3.;
    _I2MP = 3.;
    _JMP = 1.;

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

  // Operators


  // General methods
        
  void HomogMP::updateState(bool fl0, bool fl1, bool fl2)
  {
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

    Tensor3D C;
    C = (tvmet::trans(_F))*_F;
   
    const double I1=tvmet::trace(C),I2=1./2.*(I1*I1-tvmet::trace(C*C)),J=jac;

      //energy = C10(I1bar-I1MP)^2 + C01(I2bar-I2MP)^2 + D1(J-JMP)^2
    if(fl0)
      _W =   _C10*(I1*pow(J,-2./3.)-_I1MP)*(I1*pow(J,-2./3.)-_I1MP)
          +  _C01*(I2*pow(J,-4./3.)-_I2MP)*(I2*pow(J,-4./3.)-_I2MP)
          +  _D1*(J-_JMP)*(J-_JMP); 

    /*   if (_W < 0.0) { std::cout << _W << std::endl;
      std::cout << J << " " << _JMP << std::endl;
      std::cout << I1 << " " << _I1MP << std::endl;
      std::cout << I2 << " " << _I2MP << std::endl;
      assert(0);}
    */
/*      if(_W<-1e-10) {std::cout <<"Energy density is negative, something wrong wrong wrong"<<std::endl;
                std::cout <<  _W << std::endl
                          << I1*pow(J,-2./3.)-3. <<"\t"<< I2*pow(J,-4./3.)-3. <<"\t"<< (J-1.)*(J-1.) << std::endl
                          << _C10  << "\t" << _C01 << "\t" << _D1 <<std::endl;}*/
     //1stPK stress = partial W/partial F = partial W/partial I1 * partial I1/partial F + ...(I2) + ...(J)
    if(fl1)
      _P =   _C10*2*(I1*pow(J,-2./3.)-_I1MP)*pow(J,-2./3.)*(2.*_F)
          +  _C01*2*(I2*pow(J,-4./3.)-_I2MP)*pow(J,-4./3.)*(2.*_F*I1-2.*_F*C) 
          +  (   _C10*2*(I1*pow(J,-2./3.)-_I1MP)*I1*(-2./3.*pow(J,-5./3.))
               + _C01*2*(I2*pow(J,-4./3.)-_I2MP)*I2*(-4./3.*pow(J,-7./3.))
               + 2.*_D1*(J-_JMP))*J*tvmet::trans(invF);
  
    if (fl2) {
    }

    return;
  }

  double HomogMP::vonMisesStress() const{
      double jac = determinant(_F);

      Tensor3D transF(0.0);
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  transF(i,J) = _F(J,i);
	}
      }

      Tensor3D Cauchy(0.0);
      // = 1.0/jac*P*transF;

      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  Cauchy(i,J) = 1.0/jac*_P(i,J)*transF(i,J);
	}
      }

      double vonMises = 0.7071068*sqrt((Cauchy(0,0)-Cauchy(1,1))*(Cauchy(0,0)-Cauchy(1,1)) + 
				(Cauchy(1,1)-Cauchy(2,2))*(Cauchy(1,1)-Cauchy(2,2)) + 
				(Cauchy(2,2)-Cauchy(0,0))*(Cauchy(2,2)-Cauchy(0,0)) + 
				6.0*(Cauchy(0,1)*Cauchy(0,1)+Cauchy(0,2)*Cauchy(0,2)+Cauchy(1,2)*Cauchy(1,2)));
      return vonMises;
  }
    

  void HomogMP::ConsistencyTest()
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
