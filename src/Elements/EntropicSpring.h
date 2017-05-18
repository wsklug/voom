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
#include "Spring.h"
#include "VoomMath.h"

namespace voom
{

  template<int N>
  class EntropicSpring : public Spring<N> {
    
  public: 

    using Spring<N>::_nodeA;
    using Spring<N>::_nodeB;
    using Spring<N>::_k;
    using Spring<N>::_d0;
    using Spring<N>::_energy;
    using Spring<N>::_baseNodes;

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;

    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double kT, double kC,
		   int fitOrder, bool coeffsKnown, double maxForce)
                  : Spring<N> (nodeA, nodeB, kC ), _coeffsKnown(coeffsKnown), 
		    _fitOrder(fitOrder), _mult(1) , _maxForce(maxForce) { 
      const VectorND & XA = _nodeA->position();
      const VectorND & XB = _nodeB->position();
      _Lp = _k/kT; //_k = kC
      _d0 = norm2(XB-XA)/_Lp;
      _xarc = _d0*(1.0+(_d0/6.0));
      
      _fitCoeffs=0;
      
      if(_coeffsKnown) {
        getCoeffs();
      }
      
      getMaxLen();
    }

    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double kT, double kC, 
		   int fitOrder, int mult, double Larc, bool coeffsKnown, double maxForce)
                  : Spring<N> (nodeA, nodeB, kC ), _coeffsKnown(coeffsKnown), 
		    _fitOrder(fitOrder), _mult(mult), _maxForce(maxForce) {
      _Lp = _k/kT;
      _xarc = Larc/_Lp;
      _d0 = 3.0*(sqrt(1.0+(2.0*_xarc)/3.0)-1.0);
      
      _fitCoeffs=0;
      
      if(_coeffsKnown) {
        getCoeffs();
      }

      getMaxLen();
    }

    virtual void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double Lsep  = _mult*norm2(xB-xA);
      if(Lsep > _maxLen) Lsep = _maxLen;

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
    
    void getCoeffs() { // return coefficients for either a 1st, 3rd, or 5th order polynomial fit (+ divergent term) //
      _fitCoeffs = new double[1+_fitOrder];
      _fitCoeffs[0]=(3.0*_k*sqr(1.0+(_d0/6.0)))/(2.0*_Lp*_Lp);
      _fitCoeffs[1]=(_k/(_Lp*_Lp))*(3.0*(24.0 + 8.0*_d0 - _d0*_d0))/(2.0*pow(_d0,4));
      if(_fitOrder>1) {
        _fitCoeffs[2]=(_k/(2.0*_Lp*_Lp))*(18.0*(96.0 - 48.0*_d0 - 29.0*_d0*_d0))/(7.0*pow(_d0,6));
        _fitCoeffs[3]=(_k/(6.0*_Lp*_Lp))*(324.0*(408.0 - 264.0*_d0 + 3.0*_d0*_d0 + 40.0*pow(_d0,3)))/(49.0*pow(_d0,8));
      }
      if(_fitOrder>3) {
        _fitCoeffs[4]=(_k/(24.0*_Lp*_Lp))*(2592.0*(80064.0 - 18912.0*_d0 + 14799.0*_d0*_d0 + 1800.0*pow(_d0,3) - 1775.0*pow(_d0,4)))/
        (3773.0*pow(_d0,10));
        _fitCoeffs[5]=(_k/(120.0*_Lp*_Lp))*(1166400.0*(86136.0 + 27912.0*_d0 + 21091.0*_d0*_d0 + 240.0*pow(_d0,3) - 415.0*pow(_d0,4) + 280.0*pow(_d0,5)))/
        (49049.0*pow(_d0,12));
      }
    }
    
    void destroyCoeffs() {
      delete [] _fitCoeffs;
      _fitCoeffs=0;
    }
    
    double getEnergy(double x) { // return the energy for a given filament length //
      double enerMag=0.0;
      if(!_coeffsKnown) {
        getCoeffs();
      }
      
      enerMag = _Lp*_fitCoeffs[0]*(((x-_d0)/(_xarc-x))+log((_xarc-x)/(_xarc-_d0)));
      for(int j=1;j<=_fitOrder;j++) {
        enerMag += _Lp*_fitCoeffs[j]*pow(x-_d0,j+1)/(j+1);
      }
      
      if(!_coeffsKnown) {
        destroyCoeffs();
      }
      
      return enerMag;
    }
    
    double getForceMag(double x) { // return the force for a given filament length //
      double forceMag=0.0;
      if(!_coeffsKnown) {
        getCoeffs();
      }
      
      forceMag = (_fitCoeffs[0]*(x-_d0))/sqr(_xarc-x);
      for(int j=1;j<=_fitOrder;j++) {
        forceMag += _fitCoeffs[j]*pow(x-_d0,j);
      }
//       if(forceMag > _maxForce) {
// 	forceMag = _maxForce;
//       }
//       else if(forceMag < -_maxForce) {
// 	forceMag = -_maxForce;
//       }

      if(!_coeffsKnown) {
        destroyCoeffs();
      }
      
      return forceMag;
    }

    void getMaxLen() {
      double lenLow = _d0*_Lp;
      double lenHigh = _xarc*_Lp;
      double curF = getForceMag((lenLow+lenHigh)/(2.0*_Lp));
      while(fabs(curF-_maxForce) > 1.0e-3) {
	if(curF > _maxForce) lenHigh = (lenHigh+lenLow)/2.0;
	else lenLow = (lenHigh+lenLow)/2.0;
	curF = getForceMag((lenLow+lenHigh)/(2.0*_Lp));
      }
      _maxLen = (lenLow+lenHigh)/2.0;
    }

  protected:
    
    int _fitOrder;
    int _mult;
    
    //Node_t * _nodeA;
    //Node_t * _nodeB;
    double * _fitCoeffs;
    double _Lp;
    //double _kC;
    double _xarc;
    //double _x0;
    double _maxForce;
    double _maxLen;
    bool _coeffsKnown;
    
  };
};

#endif // __EntropicSpring_h__
