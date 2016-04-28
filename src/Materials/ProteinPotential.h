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
  \file ProteinPotential.h
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
    ProteinNode(DeformationNode<3> *Host, double Length): _host(Host), _length(Length) {};

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
      double DeltaZ = fabs(a(2) - b(2));
      double DeltaZperiodic = fabs(_length - DeltaZ); // 1) if _length < 0 then no periodic BC; 2) assume peridic BC are in Z
      if (DeltaZperiodic < DeltaZ) 
	{ a(2) = 0.0; b(2) = DeltaZperiodic; };
      return tvmet::norm2(a - b);
    }

    DeformationNode<3> * getHost() { return _host; };
    void setHost(DeformationNode<3> * NewHost) { _host = NewHost; };
    
  protected:
    DeformationNode<3> * _host;
    double _length;
  };
  

  
  class ProteinPotential
  {
  public:
    virtual double computeEnergy(ProteinNode * A,  ProteinNode *B) = 0;
    virtual double computeForce(ProteinNode * A,  ProteinNode *B) = 0;
    virtual double computedWdEqPar(ProteinNode * A,  ProteinNode *B) = 0;
    virtual double computeddWddEqPar(ProteinNode * A,  ProteinNode *B) = 0;
    virtual void setEquilibriumParam(double a) = 0;
    virtual double getEquilibriumParam() = 0;
    virtual void setEquilibriumR(double R) = 0;
    virtual double getEquilibriumR() = 0;
  }; 
  
}
#endif // _PROTEINPOTENTIAL_H_
