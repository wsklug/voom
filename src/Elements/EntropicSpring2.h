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


    EntropicSpring(Node_t * nodeA, Node_t * nodeB, double kap, double Lp, double L0, double mu) : Spring<N> (nodeA,nodeB,kap), _mu(mu) {
      
      _Lp = Lp;
      _L0 = L0;
      _Larc = _L0*(1.0 + (_L0/(6.0*_Lp)));
      VectorND sep;
      sep = nodeA->position() - nodeB->position();
      double segd = norm2(sep);
      _mult = _L0/segd;
      
      getCoeffs();

      getMaxLen();
      getMinLen();

    }



    virtual void compute(bool f0, bool f1, bool f2) {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double Lsep  = norm2(xB-xA);
      double scaleSep = _mult*Lsep;
      
      if(f0) {
        _energy = getEnergy(scaleSep);
      }
      
      if(f1) {
        double fmag = getForceMag(scaleSep);
	for(int i=0; i<N; i++) {	  
	  double f = fmag*(xA(i)-xB(i))/Lsep;
	  _nodeA->addForce(i, f);
	  _nodeB->addForce(i,-f);
	}
      }

      return;
    }

    double stiffnessChange() const {
      return (stiffness() - stiffness(_L0))/stiffness(_L0);

    }

    double stiffness() const {
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double Lsep  = norm2(xB-xA);
      double scaleSep = _mult*Lsep;

      return stiffness(scaleSep);
    }

    double stiffness(double L) const {
      double stiff=0.0;

      double x = L/_Lp;
      double x0 = _L0/_Lp;
      double xarc = _Larc/_Lp;
      double xmax = _maxLen/_Lp;
      double xmin = _minLen/_Lp;
      
      if(x > xmax || x < xmin) stiff = _mu/_Larc;
      else {
	stiff = _fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-x) + _fitCoeffs[2]/sqr(x) + 2.0*(x-x0)*(_fitCoeffs[1]/pow(xarc-x,3) - _fitCoeffs[2]/pow(x,3));
      }
	
      return stiff;
    }

    double initialStiff() const {
      return stiffness(_L0);
    }

    std::vector<double> check(double L) {
      std::vector<double> efandk;
      efandk.push_back(getEnergy(L));
      efandk.push_back(getForceMag(L));
      efandk.push_back(getStiffness(L));

      return efandk;
      
    }

    void reportSpringProperties() {
      std::cout << "Entropic spring with: " << std::endl
		<< "L0 = " << _L0 << std::endl
		<< "Larc = " << _Larc << std::endl
		<< "Lmin = " << _minLen << std::endl
		<< "Lmax = " << _maxLen << std::endl;
    }

    bool checkConsistency() {
      
      bool passed = true;

      bool verbose=false;
      //bool verbose=true;

      double h = _L0*(1.0e-9);
      
      const VectorND & xA = _nodeA->point();
      const VectorND & xB = _nodeB->point();
      double Lsep  = norm2(xB-xA);
      double scaleSep = _mult*Lsep;
      
      
      double Elow = getEnergy(scaleSep);
      double Ehigh = getEnergy(scaleSep+h);

      double measForce = (Ehigh-Elow)/h;
      
      
      double realForce = getForceMag(scaleSep);
      
      double forceDiff = realForce - measForce;
      double forceAvg = .5*(realForce+measForce);
      if(forceDiff/forceAvg > 1.0e-5) {
	std::cout << "Force consistency check:" << std::endl
		  << "Output Force = " << realForce << "," << std::endl
		  << "Energy Diff'd Force = " << measForce << std::endl;
	
	passed = false;
      }
      

      
      double Flow = getForceMag(scaleSep);
      double Fhigh = getForceMag(scaleSep+h);
      
      double measStiff = (Fhigh-Flow)/h;
      
      double realStiff = stiffness();
      
      double stiffDiff = realStiff - measStiff;
      double stiffAvg = .5*(realStiff+measStiff);
      if(stiffDiff/stiffAvg > 1.0e-5) {
	std::cout << "Stiff consistency check:" << std::endl
		  << "Output Stiffness = " << realStiff << "," << std::endl
		  << "Force Diff'd Stiffness = " << measStiff << std::endl;
	
	passed = false;
      }
      
      
      
      return passed;
    }
    
  private:
    
    void getCoeffs() {
	_fitCoeffs = new double[3];
	_fitCoeffs[0] = (_k*_Lp/(pow(_L0,4)))*(36.0 - (18.0+sqr(M_PI))*(_L0/_Lp) - 1.5*sqr(_L0/_Lp));
	_fitCoeffs[1] = (1.5*_k/pow(_Lp,3))*sqr(1.0+(_L0/(6.0*_Lp)));
	_fitCoeffs[2] = _k*sqr(M_PI)/(sqr(_Lp)*_L0);
	     
    }
    
    double getEnergy(double L) { // return the energy for a given filament length //
      double enerMag=0.0;

      // simple form //
      double x = L/_Lp;
      double x0 = _L0/_Lp;
      double xarc = _Larc/_Lp;
      if(x >= x0) {
	if(x < _maxLen/_Lp) {
	  enerMag = (_Lp*_fitCoeffs[0]/2.0)*sqr(x-x0) + _Lp*_fitCoeffs[1]*(((x-x0)/(xarc-x)) + log((xarc-x)/(xarc-x0))) + _Lp*_fitCoeffs[2]*(log(x/x0)-(x-x0)/x);
	}
	else {
	  double xmax = _maxLen/_Lp;
	  enerMag = (_Lp*_fitCoeffs[0]/2.0)*sqr(xmax-x0) + _Lp*_fitCoeffs[1]*(((xmax-x0)/(xarc-xmax)) + log((xarc-xmax)/(xarc-x0))) + _Lp*_fitCoeffs[2]*(log(xmax/x0)-(xmax-x0)/xmax);
	  enerMag += (_mu/(2.0*xarc))*sqr(x-xmax);
	  enerMag += _forceatMaxLen*(x-xmax);
	}
      }
      
      else {
	if(x > _minLen/_Lp) {
	  enerMag = (_Lp*_fitCoeffs[0]/2.0)*sqr(x-x0) + _Lp*_fitCoeffs[1]*(((x-x0)/(xarc-x)) + log((xarc-x)/(xarc-x0))) + _Lp*_fitCoeffs[2]*(log(x/x0)-(x-x0)/x);
	}
	else {
	  double xmin = _minLen/_Lp;
	  enerMag = (_Lp*_fitCoeffs[0]/2.0)*sqr(xmin-x0) + _Lp*_fitCoeffs[1]*(((xmin-x0)/(xarc-xmin)) + log((xarc-xmin)/(xarc-x0))) + _Lp*_fitCoeffs[2]*(log(xmin/x0)-(xmin-x0)/xmin);
	  enerMag += (_mu/(2.0*xarc))*sqr(x-xmin);
	  enerMag += _forceatMinLen*(xmin-x);
	}
      }
      
      enerMag *= _Lp;
      return enerMag;
    }
    
    double getForceMag(double L) { // return the force for a given filament length //
      double forceMag=0.0;

      double x = L/_Lp;
      double x0 = _L0/_Lp;
      double xarc = _Larc/_Lp;
      
      if(x>=x0) {
	if(x< _maxLen/_Lp) {
	  forceMag = _Lp*(x-x0)*(_fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-x) + _fitCoeffs[2]/sqr(x));
	}
	else {
	  double xmax = _maxLen/_Lp;
	  forceMag = _Lp*(xmax-x0)*(_fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-xmax) + _fitCoeffs[2]/sqr(xmax));
	  forceMag += (_mu/_Larc)*(L-_maxLen);
	}
      }
      
      else {
	if(x > _minLen/_Lp) {
	  forceMag = _Lp*(x-x0)*(_fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-x) + _fitCoeffs[2]/sqr(x));
	}
	else {
	  double xmin = _minLen/_Lp;
	  forceMag = _Lp*(xmin-x0)*(_fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-xmin) + _fitCoeffs[2]/sqr(xmin));
	  forceMag += (_mu/_Larc)*(L-_minLen);
	}
      }
      
      return forceMag;
    }

    double getStiffness(double L) {
      double stiff=0.0;

      double x = L/_Lp;
      double x0 = _L0/_Lp;
      double xarc = _Larc/_Lp;
      
      stiff = _fitCoeffs[0] + _fitCoeffs[1]/sqr(xarc-x) + _fitCoeffs[2]/sqr(x) + 2.0*(x-x0)*(_fitCoeffs[1]/pow(xarc-x,3) - _fitCoeffs[2]/pow(x,3));
	
      return stiff;
    }

    double getLength() {
      return _L0/_mult;
    }

    double getSegLength() {
      return _L0/_mult;
    }
    
    void getMaxLen() {
      double maxStiff = _mu/_Larc;
      if(getStiffness(_L0) > maxStiff) {
	_maxLen = _L0;
	_forceatMaxLen = getForceMag(_L0);
      }
      else {
	double lenLow = _L0;
	double lenHigh = _Larc;
	double curstiff = getStiffness((lenLow+lenHigh)/2.0);
	while(fabs(curstiff-maxStiff)/maxStiff > 1.0e-6 && fabs(lenHigh-lenLow)/lenLow > 1.0e-6) {
	  if(curstiff > maxStiff) lenHigh = (lenHigh+lenLow)/2.0;
	  else lenLow = (lenHigh+lenLow)/2.0;
	  curstiff = getStiffness((lenLow+lenHigh)/2.0);
	}
	_maxLen = (lenLow+lenHigh)/2.0;
	
	_forceatMaxLen = getForceMag(_maxLen);
      }

    }

    void getMinLen() {
      double maxStiff = _mu/_Larc;
      if(getStiffness(_L0) > maxStiff) {
	_minLen = _L0;
	_forceatMinLen = -getForceMag(_L0);
      }
      else {
	double lenLow = 0.0;
	double lenHigh = _L0;
	double curstiff = getStiffness((lenLow+lenHigh)/2.0);
	while(fabs(curstiff-maxStiff)/maxStiff > 1.0e-6 && fabs(lenHigh-lenLow)/lenLow > 1.0e-6) {
	  if(curstiff > maxStiff) lenLow = (lenHigh+lenLow)/2.0;
	  else lenHigh = (lenHigh+lenLow)/2.0;
	  curstiff = getStiffness((lenLow+lenHigh)/2.0);
	}
	_minLen = (lenLow+lenHigh)/2.0;
	
	_forceatMinLen = -getForceMag(_minLen);
      }
    }

  protected:
    
    double _L0;
    double _mult;
    
    //Node_t * _nodeA;
    //Node_t * _nodeB;
    double _Lp;
    //double _kC;
    double _Larc;
    //double _x0;
    double _maxLen;
    double _forceatMaxLen;
    double _minLen;
    double _forceatMinLen;
    
    double _mu;

    double* _fitCoeffs;
    
  };
};

#endif // __EntropicSpring_h__
