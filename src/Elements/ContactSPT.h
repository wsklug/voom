// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//

#if !defined(__ContactSPT_h__)
#define __ContactSPT_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Constraint.h"

namespace voom
{
/*   class Contact */
/*   { */
/*   public: */

/*     virtual void predict() = 0; */
/*     virtual void correct() = 0; */

/*   }; */

  class RigidHemisphereContactSPT : public Constraint
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef XCNode<3> xcNode;
    typedef std::vector< xcNode* > xcNodeContainer;
    typedef xcNodeContainer::iterator xcNodeIterator;
    typedef xcNodeContainer::const_iterator ConstxcNodeIterator;

    RigidHemisphereContactSPT( const xcNodeContainer & nodes, 
			    double R, Vector3D xc ) {
      _nodes = nodes;

      _R = R;
      _xc = xc;
      _FZ = 0.0;
      _active.resize(_nodes.size());
      _forces.resize(_nodes.size());
      _forces = Vector3D(0.0);
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

	  for (int i=0; i<3; i++)
	    _forces(a)[i] = _nodes[a]->getForce(i);

	  double fn =  dot(_forces(a),n);
	  if( fn > 0 ) { // non-stick contact; 
	    _forces(a) = -fn*n;
	    _nodes[a]->updateForce( _forces(a) );
	    _FZ += _forces(a)(2);
	  }
	} else {
	  // Node is not in contact; contact force is zero
	  _forces(a) = 0.0, 0.0, 0.0;
	}
      }
      return;
    }

    double FZ() const {return _FZ;}

  private:

    double _R;
    double _FZ;
    Vector3D _xc;
    xcNodeContainer _nodes;

    std::vector<bool> _active;
    blitz::Array<Vector3D,1> _forces;

  };

} // namespace voom

#endif // __ContactSPT_h__
