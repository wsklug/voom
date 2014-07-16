// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Entropic_Spring_h__)
#define __Entropic_Spring_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"

namespace voom
{

  template<int N>
  class EntropicSpring : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double kT, double kC, bool coeffsKnown) 
      : _nodeA(nodeA), _nodeB(nodeB), _kC(kC), _coeffsKnown(coeffsKnown) { 
      const VectorND & XA = _nodeA->position();
      const VectorND & XB = _nodeB->position();
      _Lp = _kC/kT;
      _x0 = norm2(XB-XA)/_Lp;
      _xarc = _x0*(1.0+(_x0/6.0));

      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);
      
      _fitCoeffsStretch=0;
      _fitCoeffsCompress=0;
      
      if(_coeffsKnown) {
        getCoeffs();
      }
    }

    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double kT, double kC, double Larc, bool coeffsKnown) 
      : _nodeA(nodeA), _nodeB(nodeB), _kC(kC), _coeffsKnown(coeffsKnown) {
      _Lp = _kC/kT;
      _xarc=Larc/_Lp;
      _x0 = 3.0*(sqrt(1.0+(2.0*_xarc)/3.0)-1.0);
      
      _baseNodes.push_back(_nodeA);
      _baseNodes.push_back(_nodeB);
      
      _fitCoeffsStretch=0;
      _fitCoeffsCompress=0;
      
      if(_coeffsKnown) {
        getCoeffs();
      }
    }

    void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double Lsep  = norm2(xB-xA);

      if(f0) {
        _energy = getEnergy(Lsep/_Lp);
      }
      
      if(f1) {
        double fmag = getForceMag(Lsep/_Lp);
	for(int i=0; i<N; i++) {	  
	  double f = fmag*(xA(i)-xB(i))/Lsep;
	  _nodeA->addForce(i, f);
	  _nodeB->addForce(i,-f);
	}
      }

      return;
    }
    
  private:
    
    void getCoeffs() { // compute coefficients for the tension fit //
      _fitCoeffsStretch = new double[5];
      _fitCoeffsStretch[0] = (5.0*(3.0 + _x0))/6.0;
      _fitCoeffsStretch[1] = (-5.0*(48.0 + 32.0*_x0 + 3.0*sqr(_x0)))/(28.*sqr(_x0));
      _fitCoeffsStretch[2] = (60.0*(6.0 + 20.0*_x0 + 9.0*sqr(_x0) + pow(_x0,3)))/(49.0*pow(_x0,4));
      _fitCoeffsStretch[3] = (-15.0*(-9360.0 - 5088.0*_x0 + 6696.0*sqr(_x0) + 3336.0*pow(_x0,3) + 355.0*pow(_x0,4)))/(3773.0*pow(_x0,6));
      _fitCoeffsStretch[4] = (90.0*(87120.0 + 23520.0*_x0 - 11880.0*sqr(_x0) + 15384.0*pow(_x0,3) + 7985.0*pow(_x0,4) + 840.0*pow(_x0,5)))/(49049.0*pow(_x0,8));

      _fitCoeffsCompress = new double[7];
      double num,denom;
      num = 6.0*(-7.024625680045251e7 + 1.9676057012259707e9*_x0 + 3.071924968929459e9*sqr(_x0) + 1.9963779860727286e9*pow(_x0,3) + 6.699136447869978e8*pow(_x0,4) + 1.1193998982064294e8*pow(_x0,5) + 3.356503462621206e6*pow(_x0,6) - 2.217755234616453e6*pow(_x0,7) - 4.733105756466191e5*pow(_x0,8) - 6.108788006281282e4*pow(_x0,9) - 5.902924175140915e3*pow(_x0,10) + 1.8767122267691776e2*pow(_x0,11) + 1.3307833534637132e2*pow(_x0,12) + 9.834746536842325*pow(_x0,13));
      denom = 9.869604401089358*(-1.8251732957814954e6*sqr(_x0) + 7.51473421768199e7*pow(_x0,3) + 1.0656918223229192e8*pow(_x0,4) + 6.142336652702212e7*pow(_x0,5) + 1.6351884152197968e7*pow(_x0,6) + 1.441804763010045e6*pow(_x0,7) - 2.5780567139182486e5*pow(_x0,8) - 8.014949652568558e4*pow(_x0,9) - 8.777074156595784e3*pow(_x0,10) - 7.168258432328461e2*pow(_x0,11) - 2.7376998670024e1*pow(_x0,12) + 1.253966942099414e1*pow(_x0,13) + 1.485678760079027*pow(_x0,14));
      _fitCoeffsCompress[0] = num/denom;
      
      num = 300.0*(-6.938443333039173e9 - 3.705826018499381e10*_x0 - 3.114415862574134e10*sqr(_x0) - 1.0971261828862179e10*pow(_x0,3) - 1.5422331877300625e9*pow(_x0,4) + 1.3060952005380238e8*pow(_x0,5) + 8.589935343354148e7*pow(_x0,6) + 1.5013020374380086e7*pow(_x0,7) + 1.5694515522264782e6*pow(_x0,8) + 7.899956357092459e4*pow(_x0,9) - 2.0798621899677342e4*pow(_x0,10) - 6.4600146846609015e3*pow(_x0,11) - 6.943452372534567e2*pow(_x0,12) - 4.6500820626325734e1*pow(_x0,13) - 4.853915023740076*pow(_x0,14));
      denom = 69.0872308076255*(-1.8251732957814954e6*pow(_x0,4) + 7.51473421768199e7*pow(_x0,5) + 1.0656918223229192e8*pow(_x0,6) + 6.142336652702212e7*pow(_x0,7) + 1.6351884152197968e7*pow(_x0,8) + 1.441804763010045e6*pow(_x0,9) - 2.5780567139182486e5*pow(_x0,10) - 8.014949652568558e4*pow(_x0,11) - 8.777074156595784e3*pow(_x0,12) - 7.168258432328461e2*pow(_x0,13) - 2.7376998670024e1*pow(_x0,14) + 1.253966942099414e1*pow(_x0,15) + 1.485678760079027*pow(_x0,16));
      _fitCoeffsCompress[1]= num/denom;
      
      num = -60000.0*(-1.9643576073266113e10 - 2.0181339494562787e10*_x0 + 1.1758673925069587e10*sqr(_x0) + 1.6996163452635805e10*pow(_x0,3) + 6.88371556154344e9*pow(_x0,4) + 1.1739104801727167e9*pow(_x0,5) - 8.293061937847634e6*pow(_x0,6) - 4.132482610219992e7*pow(_x0,7) - 8.715256128559645e6*pow(_x0,8) - 1.1125002891842553e6*pow(_x0,9) - 8.424607580996123e4*pow(_x0,10) + 7.219259512918169e3*pow(_x0,11) + 3.329768050716951e3*pow(_x0,12) + 4.2762737896893343e2*pow(_x0,13) + 3.866143365972014e1*pow(_x0,14) + 3.53900973882659*pow(_x0,15));
      denom = 6286.938003493921*(-1.8251732957814954e6*pow(_x0,6) + 7.51473421768199e7*pow(_x0,7) + 1.0656918223229192e8*pow(_x0,8) + 6.142336652702212e7*pow(_x0,9) + 1.6351884152197968e7*pow(_x0,10) + 1.441804763010045e6*pow(_x0,11) - 2.5780567139182486e5*pow(_x0,12) - 8.014949652568558e4*pow(_x0,13) - 8.777074156595784e3*pow(_x0,14) - 7.168258432328461e2*pow(_x0,15) - 2.7376998670024e1*pow(_x0,16) + 1.253966942099414e1*pow(_x0,17) + 1.485678760079027*pow(_x0,18));
      _fitCoeffsCompress[2]= num/denom;
      
      num = .18*(-1.447834368416961e8 + 2.869150577635491e9*_x0 + 3.7071750124717254e9*sqr(_x0) + 1.8311493043770276e9*pow(_x0,3) + 3.9690395242166144e8*pow(_x0,4) + 2.8888350086110118e7*pow(_x0,5) + 4.919048375348429e4*pow(_x0,6) + 9.237640440035103e5*pow(_x0,7) + 1.9957915545302272e5*pow(_x0,8) - 2.178446527748445e4*pow(_x0,9) - 6.444500040669365e3*pow(_x0,10) + 4.608802581018795e2*pow(_x0,11) + 1.6243409187000706e2*pow(_x0,12) + 8.127150250593591*pow(_x0,13));
     denom = -1.8251732957814954e6*sqr(_x0) + 7.51473421768199e7*pow(_x0,3) + 1.0656918223229192e8*pow(_x0,4) + 6.142336652702212e7*pow(_x0,5) + 1.6351884152197968e7*pow(_x0,6) + 1.441804763010045e6*pow(_x0,7) - 2.5780567139182486e5*pow(_x0,8) - 8.014949652568558e4*pow(_x0,9) - 8.777074156595784e3*pow(_x0,10) - 7.168258432328461e2*pow(_x0,11) - 2.7376998670024e1*pow(_x0,12) + 1.253966942099414e1*pow(_x0,13) + 1.485678760079027*pow(_x0,14); 
      _fitCoeffsCompress[3]= num/denom;
      
      num = 90.0*(-2.3137915145591983e9 - 1.3328116843743644e10*_x0 - 1.1959156939076097e10*sqr(_x0) - 4.775196632972044e9*pow(_x0,3) - 9.448756259366635e8*pow(_x0,4) - 4.801595850190331e7*pow(_x0,5) + 2.029833546391032e7*pow(_x0,6) + 5.654228725521835e6*pow(_x0,7) + 8.12489558749545e5*pow(_x0,8) + 7.416480412444726e4*pow(_x0,9) - 4.03768381086847e2*pow(_x0,10) - 1.740498442339042e3*pow(_x0,11) - 3.081378914670947e2*pow(_x0,12) - 3.1592133947775828e1*pow(_x0,13) - 2.461878741085571*pow(_x0,14));
      denom = 7.0*(-1.8251732957814954e6*pow(_x0,4) + 7.51473421768199e7*pow(_x0,5) + 1.0656918223229192e8*pow(_x0,6) + 6.142336652702212e7*pow(_x0,7) + 1.6351884152197968e7*pow(_x0,8) + 1.441804763010045e6*pow(_x0,9) - 2.5780567139182486e5*pow(_x0,10) - 8.014949652568558e4*pow(_x0,11) - 8.777074156595784e3*pow(_x0,12) - 7.168258432328461e2*pow(_x0,13) - 2.7376998670024e1*pow(_x0,14) + 1.253966942099414e1*pow(_x0,15) + 1.485678760079027*pow(_x0,16));
      _fitCoeffsCompress[4]= num/denom;
      
      num = -18000.0*(-1.6347722229532312e10 - 6.224679475010416e10*_x0 - 5.836288046298794e10*sqr(_x0) - 2.60691112716778e10*pow(_x0,3) - 6.212432526051057e9*pow(_x0,4) - 6.037933136346259e8*pow(_x0,5) + 8.329589603974095e7*pow(_x0,6) + 3.904359820313761e7*pow(_x0,7) + 7.093766859850867e6*pow(_x0,8) + 8.005704498943686e5*pow(_x0,9) + 3.29975863733306e4*pow(_x0,10) - 9.673284509676145e3*pow(_x0,11) - 2.519576290525897e3*pow(_x0,12) - 3.227228267503386e2*pow(_x0,13) - 3.022785569494608e1*pow(_x0,14) - 1.8068775169722426*pow(_x0,15));
      denom = 637.0*(-1.8251732957814954e6*pow(_x0,6) + 7.51473421768199e7*pow(_x0,7) + 1.0656918223229192e8*pow(_x0,8) + 6.142336652702212e7*pow(_x0,9) + 1.6351884152197968e7*pow(_x0,10) + 1.441804763010045e6*pow(_x0,11) - 2.5780567139182486e5*pow(_x0,12) - 8.014949652568558e4*pow(_x0,13) - 8.777074156595784e3*pow(_x0,14) - 7.168258432328461e2*pow(_x0,15) - 2.7376998670024e1*pow(_x0,16) + 1.253966942099414e1*pow(_x0,17) + 1.485678760079027*pow(_x0,18));
      _fitCoeffsCompress[5]= num/denom;
      
      num = -405000.0*(2.254524275897377e11 + 4.517358266129597e11*_x0 + 4.563539554083141e11*sqr(_x0) + 2.565583936403548e11*pow(_x0,3) + 8.334544702990472e10*pow(_x0,4) + 1.437699959134924e10*pow(_x0,5) + 3.1690849740701076e8*pow(_x0,6) - 4.65918295479515e8*pow(_x0,7) - 1.2604240733885201e8*pow(_x0,8) - 1.8796730652307528e7*pow(_x0,9) - 1.6665365580007752e6*pow(_x0,10) + 2.3044811021359726e4*pow(_x0,11) + 3.965846157834984e4*pow(_x0,12) + 7.388051600755239e3*pow(_x0,13) + 8.332730249605584e2*pow(_x0,14) + 7.348009457131429e1*pow(_x0,15) + 3.9488540552972275*pow(_x0,16));
      denom = 49049.0*(-1.8251732957814954e6*pow(_x0,8) + 7.51473421768199e7*pow(_x0,9) + 1.0656918223229192e8*pow(_x0,10) + 6.142336652702212e7*pow(_x0,11) + 1.6351884152197968e7*pow(_x0,12) + 1.441804763010045e6*pow(_x0,13) - 2.5780567139182486e5*pow(_x0,14) - 8.014949652568558e4*pow(_x0,15) - 8.777074156595784e3*pow(_x0,16) - 7.168258432328461e2*pow(_x0,17) - 2.7376998670024e1*pow(_x0,18) + 1.253966942099414e1*pow(_x0,19) + 1.485678760079027*pow(_x0,20));
      _fitCoeffsCompress[6]= num/denom;
    }
    
    void destroyCoeffs() {
      delete [] _fitCoeffsStretch;
      delete [] _fitCoeffsCompress;
      _fitCoeffsStretch=0;
      _fitCoeffsCompress=0;
    }
    
    double getEnergy(double x) { // return the energy for a given filament length //
      double enerMag = 0.0;
      double disp = x-_x0;
      if(!_coeffsKnown) {
        getCoeffs();
      }
      
      if(disp>=0.0) {
        double logResult = log(1.0 - (6.0*disp)/sqr(_x0));
        double curMult = (6.0*disp)/(-6.0*disp + sqr(_x0)) + logResult;
        enerMag = curMult*_fitCoeffsStretch[0];
        curMult += 2.0*disp + (6.0*(-1.0 + disp)*disp)/(-6.0*disp + sqr(_x0)) + ((-3.0 + sqr(_x0))*logResult)/3.0;
        enerMag += curMult*_fitCoeffsStretch[1];
        curMult += ((6.0*d*(6.0*(-2.0 + d)*d + (4.0 + 3.0*d)*sqr(_x0) - pow(_x0,4)))/(6.0*d - sqr(_x0)) + sqr(_x0)*(-4.0 + sqr(_x0))*logResult)/12.0;
        enerMag += curMult*_fitCoeffsStretch[2];
        curMult += ((6.0*disp*(18.0*sqr(disp)*(-3.0 + 2.0*disp) + 3.0*disp*(-9.0 + 4.0*disp)*sqr(_x0) + 3.0*(3.0 + 2.0*disp)*pow(_x0,4) - 2.0*pow(_x0,6)))/(6.0*disp - sqr(_x0)) + pow(_x0,4)*(-9.0 + 2.0*sqr(_x0))*logResult)/108.0;
        enerMag += curMult*_fitCoeffsStretch[3];
        curMult += (6.0*disp*(108.0*pow(disp,3)*(-4.0 + 3.0*disp) + 18.0*sqr(disp)*(-8.0 + 5.0*disp)*sqr(_x0) + 6.0*disp*(-12.0 + 5.0*disp)*pow(_x0,4) + 3.0*(8.0 + 5.0*disp)*pow(_x0,6) - 5.0*pow(_x0,8)) - pow(_x0,6)*(-6.0*disp + sqr(_x0))*(-24.0 + 5.0*sqr(_x0))*logResult)/(1296.0*(6.0*disp - sqr(_x0)));
        enerMag += curMult*_fitCoeffsStretch[4];
      }
      
      else {
        // return energy for compressed case; can't get exact answer, must work on //
      }
      
      if(!_coeffsKnown) {
        destroyCoeffs();
      }
      
      return (_kC/_Lp)*enerMag;
    }
    
    double getForceMag(double x) { // return the force for a given filament length //
      double forceMag = 0.0;
      double disp = x-_x0;
      if(!_coeffsKnown) {
        getCoeffs();
      }
      
      if(disp>=0.0) { // if the filament is stretched, use the stretch fitting function //
        for(int i=1;i<=5;i++) {
          forceMag += _fitCoeffsStretch[i-1]*pow(disp,i);
        }
        forceMag *= 1.0/sqr(_xarc-x);
      }
      
      else { // if the filament is compressed, use the compressed fitting function //
        double Px = 1.0;
        double Qx = 1.0;
        for(int i=1;i<=3;i++) {
          Px += _fitCoeffsCompress[i-1]*pow(disp,i);
          Qx += _fitCoeffsCompress[i+2]*pow(disp,i);
        }
        Qx += _fitCoeffsCompress[6]*pow(disp,4);
        forceMag = (9.869604401089358/sqr(x))*((Px/Qx)-1.0);
      }

      if(!_coeffsKnown) {
        destroyCoeffs();
      }
      
      return (_kC/sqr(_Lp))*forceMag;
    }
    
    Node_t * _nodeA;
    Node_t * _nodeB;
    double * _fitCoeffsStretch;
    double * _fitCoeffsCompress;
    double _Lp;
    double _kC;
    double _xarc;
    double _x0;
    bool _coeffsKnown;
    
  };
};

#endif // __EntropicSpring_h__
