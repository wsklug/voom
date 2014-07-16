// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__Contact_h__)
#define __Contact_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Constraint.h"

namespace voom
{

  class RigidHemisphereContact : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidHemisphereContact( const DefNodeContainer & nodes, 
			    double R, Vector3D xc, double friction=0.0,
			    double forceTol=1.0e-4 ) {
      _nodes = nodes;

      _R = R;
      _R0 = R;
      _xc = xc;
      _FZ = 0.0;
      _active.resize(_nodes.size());
      _forces.resize(_nodes.size());
      _forces = Vector3D(0.0);
      _friction = friction;
      if( _friction < 0.0 ) {
	_friction = -_friction;
	std::cout << "RigidHemisphereContact(): Friction coefficient should be non-negative.  Taking absolute value: " 
		  << _friction 
		  << std::endl;
      }   
      _forceTol = forceTol;
    }

    virtual void predict() {
      Vector3D x;
      for(int a=0; a<_nodes.size(); a++) {
	x = _nodes[a]->point();
	double R =  norm2(x-_xc);
 	if( R > _R ) {
	  // node is not in contact; keep it out of the active set
	  _active[a] = false;
	} else {
	  // node is penetrating hemisphere; project to surface
	  //std::cout << "Node " << a << "is penetrating AFM tip." << std::endl;
	  _active[a] = true;
	  double dR = R - _R;
	  Vector3D dx;
	  dx = dR*(x-_xc)/R;
	  x -= dx;
	  _nodes[a]->setPoint(x);
	} 
      }
      return;
    }

    virtual void correct() {
      _FZ = 0.0;
      for(int a=0; a<_nodes.size(); a++) {
	if( _active[a] ) {
	  // Node is in contact; enforce equilibrium normal to surface
	  // (also tangential if _friction=true)
	  Vector3D n;
	  n = _nodes[a]->point() - _xc;
	  n /=  norm2(n);
	  _forces(a) = _nodes[a]->force();
	  double fn =  dot(_forces(a),n);
	  if( fn > _forceTol) { // non-stick contact;
	    // tangent force 
	    double ft = norm2( _forces(a) - fn*n );
	    // check for slip
	    if( ft < _friction*fn ) { 
	      // stick case
	      _forces(a) = -_forces(a); 
	    } else { 
	      // slip case
	      _forces(a) = -fn*n;
	    }
	    _nodes[a]->updateForce( _forces(a) );
	    _FZ += _forces(a)(2);
	  } else { // node pulling away from surface
	    _forces(a) = 0.0, 0.0, 0.0;
	  }
	} else {
	  // Node is not in contact; contact force is zero
	  _forces(a) = 0.0, 0.0, 0.0;
	}
      }
      return;
    }

    double FZ() const {return _FZ;}

    double energy() const {return -_FZ*(_R-_R0);}

    double getForce(int a, int i) const {
      assert( a >= 0 && a < _forces.size() );
      assert( i >= 0 && i < 3 );
      return _forces(a)(i);
    }

    int active() const { 
      int n=0; 
      for(int i=0; i<_active.size(); i++) {
	if (_active[i]) n++;
      }
      return n; 
    }

  private:

    double _R;
    double _R0;
    double _FZ;
    Vector3D _xc;
    DefNodeContainer _nodes;

    std::vector<bool> _active;
    blitz::Array<Vector3D,1> _forces;
    double _friction;
    double _forceTol;
  };


  class RigidPlateContact : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidPlateContact( const DefNodeContainer & nodes, 
		       double Z, bool up, double friction=0.0,
		       double forceTol=1.0e-4 ) {
      _nodes = nodes;

      _Z = Z;
      _Z0 = Z;
      _up = up;
      _FZ = 0.0;
      _active.resize(_nodes.size());
      _forces.resize(_nodes.size());
      _forces = Vector3D(0.0);
      _friction = friction;
      if( _friction < 0.0 ) {
	_friction = -_friction;
	std::cout << "RigidPlateContact(): Friction coefficient should be non-negative.  Taking absolute value: " 
		  << _friction 
		  << std::endl;
      }
      _forceTol = forceTol;
    }

    virtual void predict() {
      Vector3D x;
      for(int a=0; a<_nodes.size(); a++) {
	x = _nodes[a]->point();
	double Z =  x(2);
 	if( ( _up && Z > _Z ) || ( !_up && Z < _Z ) ) {
	  // node is not in contact; keep it out of the active set
	  _active[a] = false;
	} else {
	  // node is penetrating hemisphere; project to surface
	  //std::cout << "Node " << a << "is penetrating AFM tip." << std::endl;
	  _active[a] = true;
	  _nodes[a]->setPoint(2,_Z);
	} 
      }
      return;
    }

    virtual void correct() {
      _FZ = 0.0;
      for(int a=0; a<_nodes.size(); a++) {
	if( _active[a] ) {
	  // Node is in contact; enforce equilibrium normal to surface
	  // (also tangential if _friction=true)
	  _forces(a) = _nodes[a]->force();
	  Vector3D n;
	  if(_up) n = 0.0, 0.0,  1.0;
	  else    n = 0.0, 0.0, -1.0;
	  double fn =  dot(_forces(a),n);
	  
	  if( fn > _forceTol ) { // node pushing towards surface
	    // tangent force 
	    double ft = norm2( _forces(a) - fn*n );
	    // check for slip
	    if( ft < _friction*fn ) { 
	      // stick case
	      _forces(a) = -_forces(a); 
	    } else { 
	      // slip case
	      _forces(a) = -fn*n;
	    }
	    _nodes[a]->updateForce( _forces(a) );
	    _FZ += _forces(a)(2);
	  } else { // node pulling away from surface
	    _forces(a) = 0.0, 0.0, 0.0;
	  }
	} else {
	  // Node is not in contact; contact force is zero
	  _forces(a) = 0.0, 0.0, 0.0;
	}
      }
      return;
    }

    double FZ() const {return _FZ;}

    double energy() const {return -_FZ*(_Z-_Z0);}

    void setZ(double Z) {_Z=Z;}

    double getForce(int a, int i) const {
      assert( a >= 0 && a < _forces.size() );
      assert( i >= 0 && i < 3 );
      return _forces(a)(i);
    }

    int active() const { 
      int n=0; 
      for(int i=0; i<_active.size(); i++) {
	if (_active[i]) n++;
      }
      return n; 
    }

  private:

    double _Z;
    double _Z0;
    double _FZ;

    bool _up;

    DefNodeContainer _nodes;

    std::vector<bool> _active;
    blitz::Array<Vector3D,1> _forces;
    double _friction;
    double _forceTol;
  };

} // namespace voom

#endif // __Contact_h__
