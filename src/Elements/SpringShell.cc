// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#include "VoomMath.h"
#include "SpringShell.h"

namespace voom
{

  /*  Element implementing an angle spring between two triangles for
       shell bending.  Nodal arrangement is defined as follows:

    \verbatim

   2-------1
   | (0) / |
   |   /   |
   | / (1) |
   0-------3   \endverbatim

  */

  // Constructor
  SpringShell::SpringShell( const DefNodeContainer & defNodes, 
			    double k, double theta0 )
    : _defNodes(defNodes), _k(k), _theta0(theta0) 
  {
    assert(defNodes.size() == 4);
    assert(defNodes[0] != 0);
    assert(defNodes[1] != 0);
    assert(defNodes[2] != 0);
    assert(defNodes[3] != 0);

    for(int I=0; I<4; I++) _baseNodes.push_back( defNodes[I] );

    compute(true,true,false);
  }

  // Do mechanics on element; compute energy, forces, and/or stiffness.
  void SpringShell::compute(bool f0, bool f1, bool f2) {

    const Vector3D & r0 = _defNodes[0]->point();
    const Vector3D & r1 = _defNodes[1]->point();
    const Vector3D & r2 = _defNodes[2]->point();
    const Vector3D & r3 = _defNodes[3]->point();

    // normal to face 0
    Vector3D n0;
    n0 = cross(r0,r1) + cross(r1,r2) + cross(r2,r0);
    double twoA0 = norm2(n0);
    n0/=twoA0;

    // normal to face 1
    Vector3D n1;
    n1 = cross(r0,r3) + cross(r3,r1) + cross(r1,r0);
    double twoA1 = norm2(n1);
    n1/=twoA1;

    // cosine of angle between normals
    double x = dot(n0, n1);

    // angle
    _theta = acos(x);

    // torque
    double dEdtheta = _k*(_theta-_theta0);

    double dEdx = -dEdtheta/sqrt(1-x*x);
    
    // Projection operators
    
    Tensor3D P0;
    tensorProduct(n0,n0,P0);
    P0 = -P0;
    P0(0,0) += 1.0;
    P0(1,1) += 1.0;
    P0(2,2) += 1.0;
    P0 /= twoA0;

    Tensor3D P1;
    tensorProduct(n1,n1,P1);
    P1 = -P1;
    P1(0,0) += 1.0;
    P1(1,1) += 1.0;
    P1(2,2) += 1.0;
    P1 /= twoA1;

    if( f0 ) {
      _energy = _strainEnergy = 0.5*_k*sqr( _theta-_theta0 );
    }

    if( f1 ) { 

      Vector3D f;

      Vector3D P0n1; P0n1 = P0*n1;
      Vector3D P1n0; P1n0 = P1*n0;


      // node 0
      f = dEdx * ( cross( r1-r2, P0n1 ) + cross( r3-r1, P1n0 ) );
      _defNodes[0]->updateForce(f);

      // node 1
      f = dEdx * ( cross( r2-r0, P0n1 ) + cross( r0-r3, P1n0 ) );
      _defNodes[1]->updateForce(f);

      Vector3D r01; r01 = r0-r1;
      // node 2
      f = dEdx * cross( r01, P0n1 );
      _defNodes[2]->updateForce(f);

      // node 3
      f = dEdx * cross( -r01, P1n0 );
      _defNodes[3]->updateForce(f);

    }

    return;

  }
    
} // namespace voom
