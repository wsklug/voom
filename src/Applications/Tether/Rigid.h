#include <string>
#include <iostream>
#include <vector>
#include <tvmet/Vector.h>
#include <fstream>
#include "Node.h"
#include "Element.h"

#ifdef _R
#undef _R
#endif

namespace voom {

  class RigidHemisphere : public Element
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidHemisphere( const DefNodeContainer & nodes, 
		     double k, double R, Vector3D xc ) {
      _defNodes = nodes;
      for(ConstDefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) 
	_baseNodes.push_back(*n);

      _k = k;
      _R = R;
      _xc = xc;
    }

    // E = 0.25 * k * ( (x-xc)^2 - R^2 )^2
    // f = dE/dx = k * ( (x-xc)^2 - R^2 )*(x-xc)
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidHemisphere::compute()" << std::endl;
      if(f0) _energy = 0.0;

      for(DefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
	Vector3D x;
	x = (*n)->point();
	double R = tvmet::norm2(x-_xc);
	if( R < _R ) {
	  assert( R > 0 );
	  //std::cout << "RigidHemisphere: Yes" << std::endl;
// 	  if(f0) {
// 	    double e = tvmet::dot(x-_xc,x-_xc) - _R*_R;
// 	    e *= e;
// 	    _energy += 0.25*_k*e;
// 	  }
// 	  if(f1) {
// 	    Vector3D f; 
// 	    f = _k*( tvmet::dot(x-_xc,x-_xc) - _R*_R )*( x-_xc );
// 	    for(int i=0; i<3; i++) (*n)->addForce(i,f(i));
// 	  }
	  double dR = R - _R;
	  if(f0) {
	    _energy += 0.5*_k*dR*dR;
	  }
	  if(f1) {
	    Vector3D f; 
	    f = _k*dR*(x-_xc)/R;
	    for(int i=0; i<3; i++) (*n)->addForce(i,f(i));
	  }
	}
      }
    }
  private:

    double _k;
    double _R;
    Vector3D _xc;
    DefNodeContainer _defNodes;
  };

  class RigidPlate : public Element
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidPlate( const DefNodeContainer & nodes, double k, double Z ) {
      _defNodes = nodes;
      for(ConstDefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) 
	_baseNodes.push_back(*n);

      _k = k;
      _Z = Z;
      _FZ = 0.0;
    }

    void setZ(double Z) { _Z = Z; }

    double FZ() const {return _FZ;}

    // E = 0.5 * k * ( x3-Z )^2
    // f3 = dE/dx3 = k * ( x3-Z )
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidPlate::compute()" << std::endl;

      if(f0) _energy = 0.0;

      if(f1) _FZ = 0.0;
      for(DefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
	double x3 = (*n)->getPoint(2);
	if( x3 < _Z ) {
	  //std::cout << "RigidPlate: Yes" << std::endl;
	  if(f0) {
	    _energy += 0.5*_k*(x3-_Z)*(x3-_Z);
	  }
	  if(f1) {
	    double f3 = _k*( x3-_Z );
	    _FZ += f3;
	    (*n)->addForce(2,f3);
	  }
	}
      }
    }
  private:

    double _k;
    double _Z;
    double _FZ;
    DefNodeContainer _defNodes;
  };

  class PenaltyBC : public Element
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef tvmet::Vector<bool,3> VectorBC;
    typedef DeformationNode<3> DefNode;

    PenaltyBC( DefNode * node, Vector3D x0, VectorBC bc, double k ) {

      _defNode = node;
      _baseNodes.push_back(node);

      _k = k;
      _x0 = x0;
      _bc = bc;
      _f = 0.0, 0.0, 0.0;
    }

    // E = 0.5 * k * ( x3-Z )^2
    // f3 = dE/dx3 = k * ( x3-Z )
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidBC::compute()" << std::endl;

      if(f0) _energy = 0.0;

      Vector3D dx(0.0);
      for(int i=0; i<3; i++) {
	if(_bc(i)) dx(i) = _defNode->getPoint(i) - _x0(i);
      }

      if(f0) {
	_energy += 0.5*_k*tvmet::dot(dx,dx);
      }
      if(f1) {
	_f = 0.0, 0.0, 0.0;
	for(int i=0; i<3; i++) {
	  if(_bc(i)) {
	    _f(i) = _k*dx(i);
	    _defNode->addForce(i, _f(i));
	  }
	}
      }
    }

    double force(int i) const {return _f(i);}

  private:

    double _k;
    Vector3D _x0;
    Vector3D _f;
    VectorBC _bc;
    DefNode * _defNode;
  };

} // end namespace voom
