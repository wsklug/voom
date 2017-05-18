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
#include "UnduloidGeometry.h"
#include <gsl/gsl_sf_ellint.h>

#define PI 3.141592653589793

namespace voom
{

  UnduloidGeometry::UnduloidGeometry(double a, double c) : _a(a),_c(c){
    _m = (_c*_c - _a*_a)/2;
    _n = (_c*_c + _a*_a)/2;
    _k = (_c*_c - _a*_a)/(_c*_c);
  }
  
  Vector2D UnduloidGeometry::getCurvilinearCoords(Vector3D R ){
      double x,y,z;
      x = R(0);
      y = R(1);
      z = R(2);
            
      double phi = std::atan2(y,x);
      double u;
      if(z < 0){
          u = acos((_n-(y*y))/_m);
      }
      else{
          u = 2*PI - acos((_n-(y*y))/_m);
      }
      Vector2D u_phi(u,phi);
      return u_phi;
  }
  
  Vector3D UnduloidGeometry::getCartesianCoords(Vector2D u_phi){
      double u, phi;
      u = u_phi(0);
      phi = u_phi(1);
      
      double x, y, z, y_u;
      y_u = std::sqrt(_n - _m*std::cos(u));
      x = y_u*std::cos(u);
      y = y_u*std::sin(u);
      
      double E, F, v, k; //Elliptic integrals
      
      v = (u - PI)/2;
      
      E = gsl_sf_ellint_E(v,_k,2);
      F = gsl_sf_ellint_F(v,_k,2);
      
      z = _a*F + _c*E;
      Vector3D R(x,y,z);
      return R;
  }
  
  tvmet::Vector<Vector3D,2> 
    UnduloidGeometry::getCovariantBaseVectors(double u, double phi){
     double y_prime, z_prime, k_term;
     
     y_prime = (_m*sin(u))/(2*sqrt(_n - _m*cos(u)));
     
     k_term = sqrt(1-_k*_k*cos(u/2)*cos(u/2));
     
     z_prime = (_a/2)/(k_term) + (_c/2)*(k_term);
     
     Vector3D g_u(y_prime*cos(phi),y_prime*sin(phi),z_prime);
     Vector3D g_phi(-y_prime*sin(phi),y_prime*cos(phi),0);
     
     tvmet::Vector<Vector3D,2> base_vectors(g_u,g_phi);
     return base_vectors;
  }
 
}
