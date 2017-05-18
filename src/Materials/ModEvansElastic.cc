// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.1  2005/09/28 03:50:48  ma
// Adding cytoskeleton part to SCElasttic, using Evans model
// to claculate cytoskeleton energy
//
// Revision 1.9  2005/06/18 17:05:02  klug
// Added Gaussian curvature term to energy.
//
// Revision 1.8  2005/05/23 17:38:39  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file ModEvansElastic.cc

  \brief Material class for modeling conformational change using double welled potential.

*/
#include <iostream>
#include "ModEvansElastic.h"

namespace voom 
{
  void ModEvansElastic::updateState(bool f0, bool f1, bool f2)
  {
    // first call updateState from the SCElastic to get bending
    // contributions
    SCElastic::updateState(f0, f1, f2);

    // OK, now add the stretching response
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(int alpha=0; alpha<2; alpha++){
      for(int i=0; i<3; i++) {
	for(int J=0; J<3; J++) {
	  F(i,J) += basis(alpha)(i)*refDual(alpha)(J);	  
	}
      }
    }
    //std::cout<<"F is="<<F(0,0)<<"\t"<<F(0,1)<<"\t"<<F(0,2)
    //         <<"\n"<<F(1,0)<<"\t"<<F(1,1)<<"\t"<<F(1,2)
    //         <<"\n"<<F(2,0)<<"\t"<<F(2,1)<<"\t"<<F(2,2)<<std::endl;


    Tensor3D C(0.0),Chat(0.0);
    C = trans(F)*F;
    Chat = trans(_invhatF)*trans(F)*F*_invhatF;
    //Chat = trans(_invhatF)*_invhatF;

    _trC = C(0,0)+C(1,1)+C(2,2);
    double trCSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    _J = sqrt((_trC*_trC - trCSquare)/2.0); 

    _trChat = Chat(0,0)+Chat(1,1)+Chat(2,2);
    double trChatSquare = 0.0;
    for(int i=0; i<3; i++) {
      for(int k=0; k<3; k++) {
	trChatSquare += Chat(i,k)*Chat(k,i);
      }
    }
    _Jhat = sqrt((_trChat*_trChat - trChatSquare)/2.0);
    //std::cout << "J= "<<_J <<"\nJhat= "<<_Jhat<<std::endl; 
    if(std::abs(_J-_Jhat)>1e-3) std::cout <<"ModEvansElastic::The second invariant is different w.r.t. two states"<<std::endl;
    
    Tensor3D P(0.0);
    //P = (_mu/_J)*F 
    //  + (1.0/_J)*(_kS*(_J-1) - 0.5*_mu*_trC/sqr(_J))*(_trC*F - F*C);
    P = ( -_mu* ( _trC/_J/_J*(_trChat/_J-2.0+_c2) + (_trC/_J-2.0+_c1)*_trChat/_J/_J) + _kS*(_J-1.0) )* 1/_J*(_trC*F-F*C)
        + _mu/_J*2.0*F*(_trChat/_J-2.0+_c2) + _mu/_J*(_trC/_J-2.0+_c1)*2.0*F*_invhatF*trans(_invhatF);
 
    if( f0 ) {
      // compute strain energy
      //_Ws = (_kS*sqr(_J-1.0)/2.0 + 0.5*_mu*(_trC/_J-1.0))/_J;
      _Ws = _kS*sqr(_J-1.0)/2.0 + _mu*(_trC/_J-2.0+_c1)*(_trChat/_J-2.0+_c2);
      _Ws =_Ws/_J;  //because in shells, the integration is done over the current area
      _W += _Ws;
    }
    
    if( f1 ) {
      //
      // -------------- stress resultants --------------
      for(int alpha=0; alpha<2; alpha++) 
	_n(alpha) += P*refDual(alpha)/_J;

      _n(2) = 0.0, 0.0, 0.0;
    }
    
    return;
  }
/*
  void ModEvansElastic::ConsistencyTest()
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
  }*/
}
