// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__ContactSNNI_h__)
#define __ContactSNNI_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Constraint.h"
//Modified from Contact.h class to include the effect of spheres in SNNI
namespace voom
{

  class RigidHemisphereContactSNNI : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidHemisphereContactSNNI( const DefNodeContainer & nodes,
     			    const std::vector<double> & node_radius,
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
      _node_radius = node_radius;
      if( _friction < 0.0 ) {
	_friction = -_friction;
	std::cout << "RigidHemisphereContactSNNI(): Friction coefficient should be non-negative.  Taking absolute value: " 
		  << _friction 
		  << std::endl;
      }
      if( _friction > 1.0 )	std::cout << "RigidHemisphereContactSNNI(): Friction coefficient greater than 1, rough contact will be used." <<std::endl; 

      _forceTol = forceTol;
    }

    virtual void predict() {
      Vector3D x;
      for(int a=0; a<_nodes.size(); a++) {
	x = _nodes[a]->point();
	//change because of SNNI
	double R =  norm2(x-_xc) - _node_radius[a];
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
	  Vector3D n;
	  n = _nodes[a]->point() - _xc;
	  n /=  norm2(n);
	  _forces(a) = _nodes[a]->force();
	  double fn =  dot(_forces(a),n);

	  if( fn > _forceTol) { 
            if(_friction<=1){
            // frictional contact;
	    // tangent force 
	      double ft = norm2( _forces(a) - fn*n );
	    // check for slip
	      if( ft < _friction*fn ) { 
	      // stick case
	      _forces(a) = -_forces(a); 
	      } else { 
	      // slip case
	      _forces(a) = -fn*n-_friction*fn*( _forces(a) - fn*n )/norm2( _forces(a) - fn*n );
	      }
           }
           else _forces(a) = -_forces(a); //rough contact

	    _nodes[a]->updateForce( _forces(a) );
	    _FZ += _forces(a)(2);
	  } else { // node pulling away from surface
	    _forces(a) = 0.0, 0.0, 0.0;
	    _nodes[a]->updateForce( _forces(a) );
	  }
	} else {
	  // Node is not in contact; contact force is zero
	  _forces(a) = 0.0, 0.0, 0.0;
	  _nodes[a]->updateForce( _forces(a) );
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

    bool active_node(int i) {
     return _active[i];
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
    std::vector<double> _node_radius;
  };


  class RigidPlateContactSNNI : public Constraint
  {
  public:

    // typedefs
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidPlateContactSNNI( const DefNodeContainer & nodes,
    		       const std::vector<double> & node_radius, 
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
      _node_radius = node_radius;
      if( _friction < 0.0 ) {
	_friction = -_friction;
	std::cout << "RigidPlateContactSNNI(): Friction coefficient should be non-negative.  Taking absolute value: " 
		  << _friction 
		  << std::endl;
      }
      if( _friction > 1.0 )	std::cout << "RigidPlateContactSNNI(): Friction coefficient greater than 1, rough contact will be used." <<std::endl; 
      _forceTol = forceTol;
    }

    virtual void predict() {
      Vector3D x;
      for(int a=0; a<_nodes.size(); a++) {
	x = _nodes[a]->point();
	double Z =  x(2);
	
	//change because of SNNI
	if( ( _up && Z > _Z+_node_radius[a] ) || ( !_up && Z < _Z-_node_radius[a] ) ) {
	  // node is not in contact; keep it out of the active set
	  _active[a] = false;
	} else {
	  // node is penetrating hemisphere; project to surface
	  //std::cout << "Node " << a << "is penetrating AFM tip." << std::endl;
	  _active[a] = true;
	  //change because of SNNI
	  if(_up)
	  _nodes[a]->setPoint(2,_Z+_node_radius[a]);
	  else
	  _nodes[a]->setPoint(2,_Z-_node_radius[a]);
	} 
      }
      return;
    }

    virtual void correct() {
      _FZ = 0.0;
      for(int a=0; a<_nodes.size(); a++) {
	if( _active[a] ) {
	  // Node is in contact; enforce equilibrium
	  _forces(a) = _nodes[a]->force();
	  Vector3D n;
	  if(_up) n = 0.0, 0.0,  1.0;
	  else    n = 0.0, 0.0, -1.0;
	  double fn =  dot(_forces(a),n);
	  
	  if( fn > _forceTol) { 
            if(_friction<=1){
            // frictional contact;
	    // tangent force 
	      double ft = norm2( _forces(a) - fn*n );
	    // check for slip
	      if( ft < _friction*fn ) { 
	      // stick case
	      _forces(a) = -_forces(a); 
	      } else { 
	      // slip case
	      _forces(a) = -fn*n-_friction*fn*( _forces(a) - fn*n )/norm2( _forces(a) - fn*n );
	      }
           }
           else _forces(a) = -_forces(a); //rough contact

	    _nodes[a]->updateForce( _forces(a) );
	    _FZ += _forces(a)(2);
	  } else { // node pulling away from surface
	    _forces(a) = 0.0, 0.0, 0.0;
	    _nodes[a]->updateForce( _forces(a) );
	  }
	} else {
	  // Node is not in contact; contact force is zero
	  _forces(a) = 0.0, 0.0, 0.0;
	  _nodes[a]->updateForce( _forces(a) );
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

    bool active_node(int i) {
     return _active[i];
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
    std::vector<double> _node_radius;
  };

} // namespace voom

#endif // __ContactSNNI_h__
