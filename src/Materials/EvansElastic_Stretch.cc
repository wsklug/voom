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

  \brief Material class for modeling conformational change - stretch and stretch direction are dof

*/
#include <iostream>
#include "EvansElastic_Stretch.h"

namespace voom 
{
  void EvansElastic_Stretch::updateState(bool f0, bool f1, bool f2)
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
    
    // Inverse of shear matrix - NOT USED ANYMORE - ref to EvansElastic_Skewed
    // InvHatF(0,0) = 1.0 + _eta*sin(_theta)*cos(_theta);   InvHatF(0,1) = -_eta*cos(_theta)*cos(_theta);       InvHatF(0,2) = 0.0;
    // InvHatF(1,0) = _eta*sin(_theta)*sin(_theta);         InvHatF(1,1) = 1.0 -_eta*sin(_theta)*cos(_theta);   InvHatF(1,2) = 0.0;
    // InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;

    InvHatF(0,0) = _eta*sin(_theta)*sin(_theta) + cos(_theta)*cos(_theta)/_eta;   InvHatF(0,1) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);  
     InvHatF(0,2) = 0.0;
    InvHatF(1,0) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);   InvHatF(1,1) = _eta*cos(_theta)*cos(_theta) + sin(_theta)*sin(_theta)/_eta;  
     InvHatF(1,2) = 0.0;
    InvHatF(2,0) = 0.0;                                  InvHatF(2,1) = 0.0;                                 InvHatF(2,2) = 1.0;

  
    C = trans(F)*F;
    Chat = InvHatF*C*InvHatF;   // Only because InvHatF is symmetric if prestrech is imposed

    double trC = C(0,0) + C(1,1) + C(2,2);
    double trCSquare = 0.0;
    for(i = 0; i < 3; i++) {
      for(k = 0; k < 3; k++) {
	trCSquare += C(i,k)*C(k,i);
      }
    }

    // cout << "C in EvansElastic_Stretch.cc = " << C << endl;
      
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
      
      switch (_WcType) {
        case 1:
	  _W += 0.5*(_WcConst[0]*pow(_eta - _WcConst[2], 2.0) - _WcConst[1]*cos(6.0*_theta))/J;
	  _Wc = 0.5*(_WcConst[0]*pow(_eta - _WcConst[2], 2.0) - _WcConst[1]*cos(6.0*_theta))/J;
	break;
        case 2:
	  _W += (_WcConst[0] + _WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	  _Wc = (_WcConst[0] + _WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	break;
        case 3:
	  if (_eta < 1-_WcConst[1]) {
	    _W -= (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	    _Wc =-(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	  }
	  else if(_eta > 1+_WcConst[1]) {
	    _W += (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	    _Wc = (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	  }
	  else {
	    _W += sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	    _Wc = sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	  }
	break;
      case 4:
	  if (_eta < 1-_WcConst[1]) {
	    _W += (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	  else if(_eta > 1+_WcConst[1]) {
	    _W += ( _WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = ( _WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	  else {
	    _W += ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI))) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI))) +_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	break;
      case 5:
	  if (_eta < 1-_WcConst[1]) {
	    _W -= (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	    _Wc =-(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	  }
	  else if(_eta > 1+_WcConst[1]) {
	    _W += (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	    _Wc = (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	  }
	  else {
	    _W += sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	    _Wc = sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	  }
	break;
      case 6:
	  if (_eta < 1-_WcConst[1]) {
	    _W += (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - (1.0/_WcConst[3]), 2.0))/J;
	    _Wc = (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - (1.0/_WcConst[3]), 2.0))/J;
	  }
	  else if(_eta > 1+_WcConst[1]) {
	    _W += (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	  else {
	    _W += ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) +
		    c1 + c2*_eta + c3*pow(_eta, 2.0) + c4*pow(_eta, 3.0) )/J;
	    _Wc = ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + 
		    c1 + c2*_eta + c3*pow(_eta, 2.0) + c4*pow(_eta, 3.0) )/J;
	  }
	/*
	if(_eta > 1+_WcConst[1]) {
	  _W += (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  _Wc = (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	}
	else {
	  _W += ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) +
		  0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  _Wc = ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + 
		  0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	}
	  if (_eta < 1-_WcConst[1]) {
	    // _W += (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow((1.0/_eta) - (1.0/_WcConst[3]), 2.0))/J;
	    // _Wc = (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow((1.0/_eta) - (1.0/_WcConst[3]), 2.0))/J;
	    _W += (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = (-_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	  else if(_eta > 1+_WcConst[1]) {
	    _W += (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = (_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)) + 0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	  else {
	    // _W += ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) +
	    //	    0.25*_WcConst[2]*(pow(1.0 + _WcConst[1] - _WcConst[3], 2) + pow((1.0/(1.0 - _WcConst[1])) - (1.0/_WcConst[3]), 2.0) ) )/J;
	    // _Wc = ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + 
	    //	    0.25*_WcConst[2]*(pow(1.0 + _WcConst[1] - _WcConst[3], 2) + pow((1.0/(1.0 - _WcConst[1])) - (1.0/_WcConst[3]), 2.0) ) )/J;
	    
	    _W += ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) +
		    0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	    _Wc = ( sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + 
		    0.5*_WcConst[2]*pow(_eta - _WcConst[3], 2.0) )/J;
	  }
	*/
	break;
      }

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

      // Computing forces conjugated to the imposed stretch magnitude and direction
      Phat = _mu*C*InvHatF/sqr(J);  
      // Testing derivation
      /*
      Tensor3D Ptest(0.0), Fhat(0.0);
      invert(InvHatF, Fhat);
      Ptest = InvHatF*Fhat;
      cout << " I " << Ptest << endl << endl;
      Ptest = Fhat*InvHatF;
      cout << " I " << Ptest << endl << endl;
      cout << " ks = " << _kS << " mu = " << _mu << endl;
      cout << "J = " << J << endl << endl;
      cout << "trChat = " << trChat << endl << endl;
      cout << "F = " << F << endl << endl;
      cout << "C = " << C << endl << endl;
      cout << "Chat = " << Chat << endl << endl;
      Ptest = (trC*F - F*C);
      cout << "(trC*F - F*C) = " << Ptest << endl << endl;
      Ptest = trans(F)*(trC*F - F*C)*trans(Fhat);
      cout << "trans(F)*(trC*F - F*C)*trans(Fhat) = " << Ptest << endl << endl;
      Ptest =  trans(F)*P*trans(Fhat);
      cout << " Pc = " << Phat << " trans(F)*P*trans(Fhat) = " << Ptest << endl;
      cout << " P = " << P << endl;
      Ptest = _mu*F*InvHatF*trans(InvHatF);
      cout << "_mu*F*InvHatF*trans(InvHatF) = " << Ptest << endl;
      Ptest = trans(F)*_mu*F*InvHatF*trans(InvHatF)*trans(Fhat);
      cout << "trans(F)*_mu*F*InvHatF*trans(InvHatF)*trans(Fhat) = " << Ptest << endl;
      Ptest = trans(InvHatF)*trans(Fhat);
      cout << "trans(F)*_mu*F*InvHatF*trans(InvHatF)*trans(Fhat) = " << Ptest << endl;
      Ptest = trans(F)*_mu*F*InvHatF;
      cout << "trans(F)*_mu*F*InvHatF = " << Ptest << endl;
      Ptest =  trans(F)*P*trans(Fhat);
      cout << " Pc = " << Phat << " trans(F)*P*trans(Fhat) = " << Ptest << endl;
      */
      // cout << "Phat = " << Phat << endl; 
      // 
      // Ptest = trans(F)*P*trans(Fhat);
      // cout << "FTPe = " << Ptest << endl << endl << endl;
      // // Phat = Ptest;
      // End testing of derivation

      invhatFderEta(0,0) = sin(_theta)*sin(_theta) - cos(_theta)*cos(_theta)/(_eta*_eta);
      invhatFderEta(0,1) = (-1.0 - (1.0/(_eta*_eta)))*cos(_theta)*sin(_theta);   invhatFderEta(0,2) = 0.0;
      invhatFderEta(1,0) = (-1.0 - (1.0/(_eta*_eta)))*cos(_theta)*sin(_theta);  
      invhatFderEta(1,1) = cos(_theta)*cos(_theta) - sin(_theta)*sin(_theta)/(_eta*_eta);
      invhatFderEta(2,0) = 0.0;                       invhatFderEta(2,1) = 0.0;                       invhatFderEta(2,2) = 0.0; 

      invhatFderTheta(0,0) = 2.0*sin(_theta)*cos(_theta)*(_eta - (1.0/_eta));  
      invhatFderTheta(0,1) = (-_eta + (1.0/_eta))*(cos(_theta)*cos(_theta) - sin(_theta)*sin(_theta));    invhatFderTheta(0,2) = 0.0;
      invhatFderTheta(1,0) = (-_eta + (1.0/_eta))*(cos(_theta)*cos(_theta) - sin(_theta)*sin(_theta));
      invhatFderTheta(1,1) = 2.0*sin(_theta)*cos(_theta)*(-_eta + (1.0/_eta));                            invhatFderTheta(1,2) = 0.0;  
      invhatFderTheta(2,0) = 0.0;                     invhatFderTheta(2,1) = 0.0;                         invhatFderTheta(2,2) = 0.0; 
											   
      _stretchForce = 0.0;
      _directionForce = 0.0;
      _phiForce = 0.0;
      for (i = 0; i < 3; i++) {
        for(k = 0; k < 3; k++) {
	  _stretchForce += Phat(i,k)*invhatFderEta(i,k);
	  _directionForce += Phat(i,k)*invhatFderTheta(i,k);
        }
      }

      
      switch (_WcType) {
        case 1:
	  _directionForce += 3.0*_WcConst[1]*sin(6.0*_theta)/J;
	  _stretchForce   += _WcConst[0]*(_eta - _WcConst[2])/J;
	break;
        case 2:
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	break;
      case 3:
	if (_eta < 1-_WcConst[1]) {
	  _directionForce += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	}
	else if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	  _stretchForce   += cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI)))/J;
	}
	break;
      case 4:
	if (_eta < 1-_WcConst[1]) {
	  _directionForce += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	  _stretchForce   += 2.0*_WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	  _stretchForce   += 2.0*_WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI))/J;
	  _stretchForce   += (cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI))) + 2.0*_WcConst[2]*(_eta - _WcConst[3]))/J;
	}
	break;
      case 5:
	if (_eta < 1-_WcConst[1]) {
	  _directionForce += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	}
	else if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	}
	break;
      case 6:
	if (_eta < 1-_WcConst[1]) {
	  _directionForce += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi) )/J;
	  _phiForce       += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi) )/J;
	  _stretchForce   += _WcConst[2]*(_eta - (1.0/_WcConst[3]))/J;
	}
	else if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += _WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += (cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) +
			      c2 + 2.0*c3*_eta + 3.0*c4*pow(_eta, 2.0) )/J;
	}
	/*
	if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += _WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += (cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + _WcConst[2]*(_eta - _WcConst[3]) )/J;
	}
	if (_eta < 1-_WcConst[1]) {
	  _directionForce += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi) )/J;
	  _phiForce       += 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi) )/J;
	  // _stretchForce   -= _WcConst[2]*((1.0/_eta) - (1.0/_WcConst[3]))*pow(_eta, -2.0)/J;
	  _stretchForce   += _WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else if (_eta > 1+_WcConst[1]) {
	  _directionForce -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= 6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _stretchForce   += _WcConst[2]*(_eta - _WcConst[3])/J;
	}
	else {
	  _directionForce -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  _phiForce       -= sin((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*6.0*_WcConst[0]*sin(6.0*(_theta + 0.5*M_PI + _phi))/J;
	  // _stretchForce   += cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi)))/J;
	  _stretchForce   += (cos((M_PI/(2.0*_WcConst[1]))*(_eta-1.0))*(M_PI/(2.0*_WcConst[1]))*(_WcConst[0]*cos(6.0*(_theta + 0.5*M_PI + _phi))) + _WcConst[2]*(_eta - _WcConst[3]) )/J;
	}
	*/
	break;

      }
    }
    
    return;
  }



  //! Return the deformation gradient F
  const Tensor3D EvansElastic_Stretch::DefGradient()
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
  const Tensor3D & EvansElastic_Stretch::cauchyStress()
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

    InvHatF(0,0) = _eta*sin(_theta)*sin(_theta) + cos(_theta)*cos(_theta)/_eta;   InvHatF(0,1) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);  
     InvHatF(0,2) = 0.0;
    InvHatF(1,0) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);   InvHatF(1,1) = _eta*cos(_theta)*cos(_theta) + sin(_theta)*sin(_theta)/_eta;  
     InvHatF(1,2) = 0.0;
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



  const std::vector<double > EvansElastic_Stretch::invariants()
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

    InvHatF(0,0) = _eta*sin(_theta)*sin(_theta) + cos(_theta)*cos(_theta)/_eta;   InvHatF(0,1) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);  
     InvHatF(0,2) = 0.0;
    InvHatF(1,0) = (-_eta + (1.0/_eta))*cos(_theta)*sin(_theta);   InvHatF(1,1) = _eta*cos(_theta)*cos(_theta) + sin(_theta)*sin(_theta)/_eta;  
     InvHatF(1,2) = 0.0;
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
