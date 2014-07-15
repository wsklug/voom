// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__NematicProbTable_h__)
#define __NematicProbTable_h__

#include "VoomMath.h"

namespace voom {

  class NematicProbTable {
    
  public: 
    
    typedef std::vector< std::pair<double,double> > CumProbDist; 

    NematicProbTable() {}
    
    void setTable(const double b) {
      if(_pdfType == 0) {
	_param = b;
	int nSteps = 1024;
	double stp_size = M_PI/(2.0*(nSteps-1));
	_cumProbs.push_back(std::pair<double,double>(0.0,0.0));
	for(int i=1; i<nSteps; i++) {
	  double curAng = (2.0*i-1.0)*stp_size/2.0;
	  double curProb = exp(b*cos(2.0*curAng)/2.0);
	  curAng += stp_size/2.0;
	  double cumProb = _cumProbs[i-1].second;
	  cumProb += curProb;
	  _cumProbs.push_back(std::pair<double,double>(curAng,cumProb));
	}
	double finProb = _cumProbs[nSteps-1].second;
	for(int i=0; i<nSteps; i++) {
	  _cumProbs[i].second /= finProb;
	}
      }
      else if(_pdfType == 1) {
	_param = sqrt((b-1.0)/b);
      }
    }

    double findAngle(const double rn) {
      if(_pdfType == 0) {
	int nPts = _cumProbs.size();
	int low = 0;
	int high = nPts-1;
	bool fd = false;
	double ang;
	while(!fd) {
	  int cur = (low+high)/2;
	  if(_cumProbs[cur].second > rn) high = cur;
	  else low = cur;
	  if(high - low == 1) {
	    double cumProbdiff = _cumProbs[high].second - _cumProbs[low].second;
	    double angdiff = _cumProbs[high].first - _cumProbs[low].first;
	    double slp = cumProbdiff/angdiff;
	    ang = _cumProbs[low].first + (rn-_cumProbs[low].second)/slp;
	    fd = true;
	  }
	}
	return ang;
      }
      else if(_pdfType == 1) {
	return atan(_param*tan(M_PI*rn/2.0));
      }
    }
    
    void setPDF(int pdfType) {
      _pdfType = pdfType;
    }
    
  private:
    
    CumProbDist _cumProbs;
    int _pdfType;
    double _param;
    
  };
};

#endif // __NematicProbTable_h__
