// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include "RigidPlateAL.h"

namespace voom
{

  RigidPlateAL::RigidPlateAL( const DefNodeContainer & nodes, double k, 
			      double Z, bool up, double friction ) {
    _defNodes = nodes;
    _active.resize(nodes.size());
    for(ConstDefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      _baseNodes.push_back(*n);
    }
    _k = k;
    _Z = Z;
    _FZ = 0.0;
    _up = up;
    if(_up) {
      _normal = 0.0, 0.0, 1.0;
    } else {
      _normal = 0.0, 0.0, -1.0;
    }

    _penetration = 0.0;

    _mu = friction;
    _penetration = 0.0;
    _active.resize(nodes.size());
    _x_stick.resize(_defNodes.size());
    _x_stick = Vector3D(0.0);
    _forces.resize(_defNodes.size());
    _forces = Vector3D(0.0);
    _pressureMultipliers.resize(_defNodes.size());
    _pressureMultipliers = 0.0;
    _frictionMultipliers.resize(_defNodes.size());
    _frictionMultipliers = Vector3D(0.0);
    updateContact();
  }

  int RigidPlateAL::active() const { 
    int n=0; 
    for(int i=0; i<_active.size(); i++) {
      if (_active[i]) n++;
    }
    return n; 
  }

  void RigidPlateAL::updateContact() {

     
    // Now update surface normals and check for
    // active contacts
    for(int a=0; a<_defNodes.size(); a++) {
      const Vector3D & x = _defNodes[a]->point();
      double z = x(2) ;
	
      // save projection of current point on surface as frictional
      // ``stick'' point
      _x_stick(a) = x;
      _x_stick(a)(2) = _Z;

      // check for active contacts
      if( (_up && z < _Z) || (!_up && z > _Z) ) {	
	_active[a] = true;
      } else {
	_active[a] = false;
      }
	
      // Uzawa update: set multipliers equal to most recently computed
      // forces
      _pressureMultipliers(a) = dot(_forces(a), _normal);
      _frictionMultipliers(a) = _forces(a) 
	- _pressureMultipliers(a)*_normal;
	

    }

    return;
  }

  void RigidPlateAL::compute( bool f0, bool f1, bool f2 ) {
    //std::cout << "RigidPlateAL::compute()" << std::endl;

    if(f0) _energy = 0.0;
      
    if(f1) {
      _FZ = 0.0;
      _forces = Vector3D(0.0);
    }

    _penetration = 0.0;


    for(int a=0; a<_defNodes.size(); a++) {
	
      DefNode * nd = _defNodes[a];
      Vector3D x;
      x = nd->point();
	
      // Do normal contact first
	
      // normal gap function
      double gn = dot(x - _x_stick(a), _normal);
      _penetration = std::max(_penetration, -gn);

      // pressure multiplier and force
      double pn = _pressureMultipliers(a);
      double fn = pn + _k*gn;
      if( fn > 0 ) { // case 1: no contact 
	if(f0) {
	  _energy -= 0.5*sqr(pn)/_k; 
	}	  
	  
      } else { // case 2: contact

	if(f0) {
	  _energy += (pn+0.5*_k*gn)*gn; 
	}
	if(f1) {
	  _forces(a) += fn*_normal;
	}
      }
	
      // now frictional contact 

      // previous contact point = xc + radius times previous normal

      // compute vector from previous contact point to current point
      Vector3D dx;
      dx = x - _x_stick(a);

      // tangential gap, (I-n^t n)*dx
      Vector3D gt;
      gt = dx - dot(_normal,dx)*_normal;

      // friction multiplier and force
      const Vector3D pt = _frictionMultipliers(a);
      double norm_pt = norm2(pt);

      Vector3D ft; 
      ft = pt + _k*gt;
	
      double norm_ft = norm2(ft);

      if( pn > 0 ) { // case 1: no contact 
	if(f0) {
	  _energy -= 0.5*sqr(norm_pt)/_k; 
	}	  
	  
      } else if( norm_ft <= _mu*std::abs(pn) ) { // case 2: contact, stick

	if(f0) {
	  _energy += dot(pt+0.5*_k*gt, gt); 
	}
	if(f1) {
	  _forces(a) += ft - dot(ft, _normal)*_normal;
	}

      } else { // case 3: contact, slip
	  
	if(f0) {
	  _energy -= 
	    0.5*( sqr(norm_pt) - 2.0*_mu*std::abs(pn)*norm_ft + sqr(_mu*pn) )/_k; 
	}
	if(f1) {
	  _forces(a) += 
	    _mu*std::abs(pn)*( ft-dot(ft,_normal)*_normal )/norm_ft;
	}

      }

      if(f1) {
	const Vector3D & f = _forces(a); 
	_FZ += f(2);
	for(int i=0; i<3; i++) nd->addForce(i,f(i)); 
      }

    }

    return;
  }
  
  
} // namespace voom
