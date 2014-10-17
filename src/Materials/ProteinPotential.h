// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*! 
  \file ProteninPotential.h
*/

#ifndef _PROTEINPOTENTIAL_H_
#define _PROTEINPOTENTIAL_H_

#include "VoomMath.h"
#include "Node.h"

#include <set>
#include <map>
#include <math.h> 

using namespace std;

namespace voom {

  class ProteinNode 
  {
  public:
    ProteinNode(DeformationNode<3> *Host): _host(Host) {};

    virtual ~ProteinNode() {};
    
    DeformationNode<3>::Point getHostPosition() { return _host->point(); };
    double getDistance(ProteinNode *B) {
      DeformationNode<3>::Point a = _host->point(), b = B->getHostPosition();
      // double R = tvmet::norm2(a);
      // a /= tvmet::norm2(a);
      // b /= tvmet::norm2(b);
      // double theta = acos(dot(a,b));
      // return R*theta;
      // cout << "test " << tvmet::norm2(a-b) << " " << norm2(a-b) << endl;
      return tvmet::norm2(a-b);
    }

    DeformationNode<3> * getHost() { return _host; };
    void setHost(DeformationNode<3> * NewHost) { _host = NewHost; };
    
  protected:
    DeformationNode<3> * _host;
  };
  

  
  class ProteinPotential
  {
  public:
    virtual double computeEnergy(ProteinNode * A,  ProteinNode *B) = 0;
    virtual double computedWdEqPar(ProteinNode * A,  ProteinNode *B) = 0;
    virtual double computeddWddEqPar(ProteinNode * A,  ProteinNode *B) = 0;
    virtual void setEquilibriumParam(double a) = 0;
    virtual double getEquilibriumParam() = 0;
    virtual void setEquilibriumR(double R) = 0;
    virtual double getEquilibriumR() = 0;
  }; 
  
}
#endif // _PROTEINPOTENTIAL_H_
