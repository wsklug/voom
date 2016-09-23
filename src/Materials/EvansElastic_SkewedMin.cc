// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                         Ankush Aggarwal
//                          Luigi Perotti
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
  \file EvansElastic_SkewedMin.cc

  \brief Material class for modeling conformational change using single SKEWED welled potential - shear value is a dof

*/
#include <iostream>
#include "EvansElastic_SkewedMin.h"

namespace voom 
{
  void EvansElastic_SkewedMin::updateState(bool f0, bool f1, bool f2)
  {
    int i = 0, j = 0, alpha = 0, k = 0;
    // Stretching response - the only one considered in this class!
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    //! InvHatF = inverse of the deformation gradient of sheared equilibrium state
    Tensor3D F(0.0), C(0.0), Chat(0.0), InvHatF(0.0);
    for(alpha = 0; alpha < 2; alpha++){
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  F(i,j) += basis(alpha)(i)*refDual(alpha)(j);	  
	}
      }
    }
    // cout << "F is = " << F(0,0) << "\t" << F(0,1) << "\t" << F(0,2) << endl
    //                   << F(1,0) << "\t" << F(1,1) << "\t" << F(1,2) << endl
    //                   << F(2,0) << "\t" << F(2,1) << "\t" << F(2,2) << endl;

    InvHatF(0,0) = 1.0 + _eta*sin(_theta)*cos(_theta);   InvHatF(0,1) = -_eta*cos(_theta)*cos(_theta);       InvHatF(0,2) = 0.0;
    InvHatF(1,0) = _eta*sin(_theta)*sin(_theta);         InvHatF(1,1) = 1.0 -_eta*sin(_theta)*cos(_theta);   InvHatF(1,2) = 0.0;
    InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;

    C = trans(F)*F;
    Chat = trans(InvHatF)*trans(F)*F*InvHatF;

    double trC = C(0,0) + C(1,1) + C(2,2);
    double trCSquare = 0.0;
    for(i = 0; i < 3; i++) {
      for(k = 0; k < 3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
      
    double J = sqrt( 0.5*(trC*trC - trCSquare) ); 
    
    double trChat = Chat(0,0) + Chat(1,1) + Chat(2,2);
    
    // double Jhat = sqrt( 0.5*(trChat*trChat - trChatSquare) );
    // cout << "J= " << J << endl << "Jhat = " << Jhat << endl; 
    // if(abs(J - Jhat) > 1e-3) cout << " EvansElastic_Skewedmin: The second invariant is different w.r.t. to states " << endl;
    


    Tensor3D P(0.0), Phat(0.0), invhatFderEta(0.0), invhatFderTheta(0.0);
    //P = (_mu/_J)*F 
    //  + (1.0/_J)*(_kS*(_J-1) - 0.5*_mu*_trC/sqr(_J))*(_trC*F - F*C);
    P = ( (trC*F - F*C)*(-0.5*_mu*trChat/sqr(J) + _kS*(J-1.0) )
	  + _mu*F*InvHatF*trans(InvHatF) )/J;

    


 
    if( f0 ) {
      // compute strain energy
      //_Ws = (_kS*sqr(_J-1.0)/2.0 + 0.5*_mu*(_trC/_J-1.0))/_J;
      _W = 0.5*(_kS*sqr(J-1.0) + _mu*(trChat/J - 2.0) )/J; 
      _Ws = _W;
      // because in shells, the integration is done over the current area Ws/J
    }
    
    if( f1 ) {
      // Set _m and _nTC to zero even if they should be never used here
      _m(0) = 0.0, 0.0, 0.0;
      _m(1) = 0.0, 0.0, 0.0;
      _nTC(0) = 0.0, 0.0, 0.0;
      _nTC(1) = 0.0, 0.0, 0.0;
      _nTC(2) = 0.0, 0.0, 0.0;
      // compute stress resultants
      for(alpha = 0; alpha < 2; alpha++) 
	_n(alpha) = P*refDual(alpha)/J;

      _n(2) = 0.0, 0.0, 0.0;

      // Computing forces conjugated to the imposed shear magnitude and direction
      Phat = _mu*C*InvHatF/sqr(J);   
      invhatFderEta(0,0) = sin(_theta)*cos(_theta);   invhatFderEta(0,1) =-cos(_theta)*cos(_theta);   invhatFderEta(0,2) = 0.0;
      invhatFderEta(1,0) = sin(_theta)*sin(_theta);   invhatFderEta(1,1) =-sin(_theta)*cos(_theta);   invhatFderEta(1,2) = 0.0;  
      invhatFderEta(2,0) = 0.0;                       invhatFderEta(2,1) = 0.0;                       invhatFderEta(2,2) = 0.0; 

      invhatFderTheta(0,0) = _eta*(pow(cos(_theta),2.0)-pow(sin(_theta),2.0));   invhatFderTheta(0,1) = 2.0*_eta*cos(_theta)*sin(_theta);    invhatFderTheta(0,2) = 0.0;
      invhatFderTheta(1,0) =  2.0*_eta*cos(_theta)*sin(_theta);   invhatFderTheta(1,1) = _eta*(pow(sin(_theta),2.0)-pow(cos(_theta),2.0));   invhatFderTheta(1,2) = 0.0;  
      invhatFderTheta(2,0) = 0.0;                             invhatFderTheta(2,1) = 0.0;                                            invhatFderTheta(2,2) = 0.0; 
											   
      _shearForce = 0.0;
      _directionForce = 0.0;
      for (i = 0; i < 3; i++) {
        for(k = 0; k < 3; k++) {
	  _shearForce += Phat(i,k)*invhatFderEta(i,k);
	  _directionForce += Phat(i,k)*invhatFderTheta(i,k);
        }
      }


    }
    
    return;
  }



  //! Return the deformation gradient F
  const Tensor3D EvansElastic_SkewedMin::DefGradient()
  {
    int i = 0, j = 0, alpha = 0;
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0);
    for(alpha = 0; alpha < 2; alpha++){
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  F(i,j) += basis(alpha)(i)*refDual(alpha)(j);	  
	}
      }
    }
    
    /*
    Tensor3D Fskew(0.0), InvHatF(0.0);
    double eta = _nodeEta->point();
    InvHatF(0,0) = 1.0 + _eta*sin(_theta)*cos(_theta);   InvHatF(0,1) = -_eta*cos(_theta)*cos(_theta);       InvHatF(0,2) = 0.0;
    InvHatF(1,0) = _eta*sin(_theta)*sin(_theta);         InvHatF(1,1) = 1.0 -_eta*sin(_theta)*cos(_theta);   InvHatF(1,2) = 0.0;
    InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;
    Fskew = F*InvHatF;
    return Fskew;
    */

    return F;
  }



  //! Return the cauchy stress
  const Tensor3D & EvansElastic_SkewedMin::cauchyStress()
  {
    int alpha = 0, i = 0, j = 0, k = 0;
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0), InvHatF(0.0), C(0.0), Chat(0.0), P(0.0);
    for(alpha = 0; alpha < 2; alpha++){
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  F(i,j) += basis(alpha)(i)*refDual(alpha)(j);	  
	}
      }
    }
    C = trans(F)*F;

    InvHatF(0,0) = 1.0 + _eta*sin(_theta)*cos(_theta);   InvHatF(0,1) = -_eta*cos(_theta)*cos(_theta);       InvHatF(0,2) = 0.0;
    InvHatF(1,0) = _eta*sin(_theta)*sin(_theta);         InvHatF(1,1) = 1.0 -_eta*sin(_theta)*cos(_theta);   InvHatF(1,2) = 0.0;
    InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;
    Chat = trans(InvHatF)*trans(F)*F*InvHatF;

    double trC = C(0,0) + C(1,1) + C(2,2);
    double trCSquare = 0.0;
    for(i = 0; i < 3; i++) {
      for(k = 0; k < 3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }
    double J = sqrt( 0.5*(trC*trC - trCSquare) ); 

    double trChat = Chat(0,0) + Chat(1,1) + Chat(2,2);
    
    P = ( (trC*F - F*C)*(-0.5*_mu*trChat/sqr(J) + _kS*(J-1.0) )
          + _mu*F*InvHatF*trans(InvHatF) )/J;

    _cauchy = P*trans(F)/J;

    return _cauchy;
  }



  const std::vector<double > EvansElastic_SkewedMin::invariants()
  {
    std::vector<double > invariants(2, 0.0);
    int alpha = 0, i = 0, j = 0, k = 0;
    typedef tvmet::Vector< Vector3D, 2 > BasisVectors;
                
    const BasisVectors & basis = _deformedGeometry.a();
    const BasisVectors & dual = _deformedGeometry.aDual();

    const BasisVectors & refBasis = _referenceGeometry.a();
    const BasisVectors & refDual = _referenceGeometry.aDual();
    
    Tensor3D F(0.0), InvHatF(0.0), C(0.0), Chat(0.0), P(0.0);
    for(alpha = 0; alpha < 2; alpha++){
      for(i = 0; i < 3; i++) {
	for(j = 0; j < 3; j++) {
	  F(i,j) += basis(alpha)(i)*refDual(alpha)(j);	  
	}
      }
    }
    C = trans(F)*F;

    InvHatF(0,0) = 1.0 + _eta*sin(_theta)*cos(_theta);   InvHatF(0,1) = -_eta*cos(_theta)*cos(_theta);       InvHatF(0,2) = 0.0;
    InvHatF(1,0) = _eta*sin(_theta)*sin(_theta);         InvHatF(1,1) = 1.0 -_eta*sin(_theta)*cos(_theta);   InvHatF(1,2) = 0.0;
    InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;
    Chat = trans(InvHatF)*trans(F)*F*InvHatF;

    double trC = C(0,0) + C(1,1) + C(2,2);
    double trCSquare = 0.0;
    for(i = 0; i < 3; i++) {
      for(k = 0; k < 3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }

    double J = sqrt( 0.5*(trC*trC - trCSquare) ); 
    double trChat = Chat(0,0) + Chat(1,1) + Chat(2,2);
    
    invariants[0] = trChat/J;
    invariants[1] = J;

    return invariants;
  }

}
