// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2016 All Rights Reserved
//
//----------------------------------------------------------------------
//

/*! 
  \file UnduloidGeometry.h

  \brief ShellGeometry is a structure which computes and contains all
  of the geometric quantities relevant to calculations on an unduloid surface.

*/

#if !defined(__ShellGeometry_h__)
#define __ShellGeometry_h__

#include<blitz/array.h> 
#include<vector>
#include<iostream>
#include "voom.h"


namespace voom
{
  class UnduloidGeometry
  {
  private:
    double _a;
    double _c;
    double _m;
    double _n;
    double _k;
  public:
    UnduloidGeometry(double, double);
    Vector3D getCartesianCoords(Vector2D);
    Vector2D getCurvilinearCoords(Vector3D);
    tvmet::Vector<Vector3D,2> getCovariantBaseVectors(double u, double phi);    
  };
}; //namespace voom
#endif //_UnduloidGeometry_h__
