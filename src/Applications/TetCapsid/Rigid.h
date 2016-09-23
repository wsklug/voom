#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <tvmet/Vector.h>
#include <fstream>
#include "Node.h"
#include "Element.h"

#ifdef _R
#undef _R
#endif

namespace voom {

  class RigidSurface : public Element
  {
  public: 
    virtual double FZ() const = 0;
    virtual double penetration() const {return _penetration;}
  protected:
    double _penetration;
  };


  class RigidHemisphere : public RigidSurface
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
      _FZ = 0.0;
      _penetration = 0.0;
    }

    double FZ() const {return _FZ;}

    // E = 0.25 * k * ( (x-xc)^2 - R^2 )^2
    // f = dE/dx = k * ( (x-xc)^2 - R^2 )*(x-xc)
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidHemisphere::compute()" << std::endl;
      if(f0) _energy = 0.0;
      if(f1) _FZ = 0.0;
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
            _FZ += f(2);
 
	  }
	}
      }
    }
  private:

    double _k;
    double _R;
    double _FZ;
    Vector3D _xc;
    DefNodeContainer _defNodes;

  };

  class RigidPlate : public RigidSurface
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidPlate( const DefNodeContainer & nodes, double k, double Z, bool up=true ) {
      _defNodes = nodes;
      for(ConstDefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) 
	_baseNodes.push_back(*n);

      _k = k;
      _Z = Z;
      _FZ = 0.0;
      _up = up;
      _penetration = 0.0;
    }

    void setPenaltyCoefficient( double k ) { _k = k; }

    double penaltyCoefficient() const {return _k;}

    void setZ(double Z) { _Z = Z; }

    double FZ() const {return _FZ;}

    double penetration() const {return _penetration;}

    // E = 0.5 * k * ( x3-Z )^2
    // f3 = dE/dx3 = k * ( x3-Z )
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidPlate::compute()" << std::endl;

      if(f0) _energy = 0.0;

      if(f1) _FZ = 0.0;

      _penetration = 0.0;
      for(DefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
	double x3 = ( (*n)->getPoint(2) );
 	if( (_up && x3 < _Z) || (!_up && x3 > _Z) ) {
	  //std::cout << "RigidPlate: Yes" << std::endl;
	  double dx3=x3-_Z;
	  _penetration = std::max( _penetration, std::abs(dx3) );
	  if(f0) {
	    _energy += 0.5*_k*dx3*dx3;
	  }
	  if(f1) {
	    double f3 = _k*dx3;
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
    bool _up;

    double _penetration;
  };

  class RigidHemisphereAL : public RigidSurface
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;

    RigidHemisphereAL( const DefNodeContainer & nodes, 
		     double k, double R, Vector3D xc ) {
      _defNodes = nodes;
      for(ConstDefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) 
	_baseNodes.push_back(*n);

      _k = k;
      _R = R;
      _xc = xc;
      _FZ = 0.0;
      _penetration = 0.0;
      _active.resize(nodes.size());
      _forces.resize(_defNodes.size());
      _forces = Vector3D(0.0);
      computeMultipliers();
    }

    double penaltyCoefficient() const {return _k;}

    void setPenaltyCoefficient( double k ) { _k = k; }

    void computeMultipliers() {
      Vector3D x;
      for(int a=0; a<_defNodes.size(); a++) {
	x = _defNodes[a]->point();
	double R = tvmet::norm2(x-_xc);
 	if( R < _R ) {
	  _active[a] = true;
	  double dR = R - _R;
	  Vector3D dx;
	  dx = dR*(x-_xc)/R;
	  _forces(a) += _k*dx;
	} else {
	  _active[a] = false;
	  _forces(a) = Vector3D(0.0);
	}
      }
      return;
    }

    double FZ() const {return _FZ;}

    double penetration() const {return _penetration;}

    // E = 0.25 * k * ( (x-xc)^2 - R^2 )^2
    // f = dE/dx = k * ( (x-xc)^2 - R^2 )*(x-xc)
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidHemisphere::compute()" << std::endl;
      if(f0) _energy = 0.0;

      if(f1) _FZ = 0.0;

      _penetration = 0.0;
 
//       for(DefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
      for(int a=0; a<_defNodes.size(); a++) {
	Vector3D x;
	x = _defNodes[a]->point();
	double R = tvmet::norm2(x-_xc);
	Vector3D f; 
	f = _forces(a);
	if( R < _R ) {
	  assert( R > 0 );
 	 _active[a] = true;
 	}
	else if( tvmet::norm2(f) == 0.0 ) _active[a] = false;
      }	
      for(int a=0; a<_defNodes.size(); a++) {
	if( _active[a] ) {
	  DefNode * n = _defNodes[a];
	  Vector3D x;
	  x = n->point();
	  double R = tvmet::norm2(x-_xc);
	  Vector3D f; 
	  f = _forces(a);
	  double dR = R - _R;
	  Vector3D dx;
	  dx = dR*(x-_xc)/R;
	  double e = tvmet::dot(f,dx);
	  if( std::abs(dR) > std::abs(_penetration) ) _penetration = dR;
	  e += 0.5*_k*dR*dR;
	  f += _k*dx;
	  if(f0) {
// 	    _energy += 0.5*_k*dR*dR;
	    _energy += e;
	  }
	  if(f1) {
// 	    Vector3D f; 
// 	    f = _k*dR*(x-_xc)/R;
            _FZ += f(2);
	    for(int i=0; i<3; i++) n->addForce(i,f(i)); 
	  }
	}
      }
    }
  private:

    double _k;
    double _R;
    double _FZ;
    Vector3D _xc;
    DefNodeContainer _defNodes;

    std::vector<bool> _active;
    blitz::Array<Vector3D,1> _forces;

  };

  class RigidPlateAL : public RigidSurface
  {
  public:

    // typedefs
    typedef tvmet::Vector<double,3> Vector3D;
    typedef DeformationNode<3> DefNode;
    typedef std::vector< DefNode* > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeIterator;
    typedef std::map< DefNode*, double > ActiveSet;
    typedef ActiveSet::iterator ActiveSetIterator;
    typedef ActiveSet::const_iterator ConstActiveSetIterator;

    RigidPlateAL( const DefNodeContainer & nodes, double k, double Z, bool up=true ) {
      _defNodes = nodes;
      _active.resize(nodes.size());
      for(ConstDefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
	_baseNodes.push_back(*n);
      }
      _k = k;
      _Z = Z;
      _FZ = 0.0;
      _up = up;
      _penetration = 0.0;
      _forces.resize(_defNodes.size());
      _forces = 0.0;
      computeMultipliers();
      compute(true,true,false);
    }

//     void detectContact() {
//       // add nodes which are penetrating and increment forces
//       for(DefNodeIterator n=_defNodes.begin(); n!=_defNodes.end(); n++) {
// 	double x3 = ( (*n)->getPoint(2) );
//  	if( (_up && x3 < _Z) || (!_up && x3 > _Z) ) {
// 	  double f3 = _k*(x3-_Z);
// 	  ActiveSetIterator a = _active.find( *n );
// 	  if( a == _active.end() ) {
// 	    // add new node
// 	    _active.insert( std::make_pair( *n, f3 ) );
// 	  } else {
// 	    // increment forces of previously penetrating nodes
// 	    a->second += f3;
// 	  }
// 	} else {
// 	  // remove nodes which are no longer penetrating	
// 	  ActiveSetIterator a = _active.find( *n );
// 	  if( a != _active.end() ) {
// 	    _active.erase( a );
// 	  }
// 	}
//       }
//       std::cout << "Detected " << _active.size() << " penetrating nodes." 
// 		<< std::endl;
//       return;
//     }

    double penaltyCoefficient() const {return _k;}

    void setPenaltyCoefficient( double k ) { _k = k; }

    void computeMultipliers() {
      for(int a=0; a<_defNodes.size(); a++) {
	double x3 = _defNodes[a]->getPoint(2) ;
 	if( (_up && x3 < _Z) || (!_up && x3 > _Z) ) {
	  _active[a] = true;
	  _forces(a) += _k*(x3-_Z);
	} else {
	  _active[a] = false;
	  _forces(a) = 0.0;
	}
      }
      return;
    }

    void setZ(double Z) { _Z = Z; }

    double FZ() const {return _FZ;}

    double penetration() const {return _penetration;}

    // E = 0.5 * k * ( x3-Z )^2
    // f3 = dE/dx3 = k * ( x3-Z )
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidPlateAL::compute()" << std::endl;

      if(f0) _energy = 0.0;
      
      if(f1) _FZ = 0.0;

      _penetration = 0.0;
      for(int a=0; a<_defNodes.size(); a++) {
	DefNode * n = _defNodes[a];
	double x3 = n->getPoint(2) ;
	double f3 = _forces(a);
 	if( (_up && x3 < _Z) || (!_up && x3 > _Z) ) _active[a] = true;
 	else if( f3 == 0.0 ) _active[a] = false;
      }
      for(int a=0; a<_defNodes.size(); a++) {
	if( _active[a] ) {
	  DefNode * n = _defNodes[a];
	  double x3 = n->getPoint(2) ;
	  double dx3 = x3 - _Z;
	  // multiplier part
	  double f3 = _forces(a);
	  double e = f3*dx3; 
// 	  if( _up ) _penetration = std::max(-dx3,_penetration);
// 	  else  _penetration = std::max(dx3,_penetration);
	  if( std::abs(dx3) > std::abs(_penetration) ) _penetration = dx3;
	  // penalty part
//	  if( /* f3 != 0.0 || */ (_up && x3 < _Z) || (!_up && x3 > _Z) ) {
	    e += 0.5*_k*dx3*dx3;
	    f3 += _k*dx3;
// 	  }

	  if(f0) {
// 	    _energy += f3*dx3 + 0.5*_k*dx3*dx3;
	    _energy += e;
	  }
	  if(f1) {
// 	    f3 += _k*dx3;
	    _FZ += f3;
	    n->addForce(2,f3);
	  }
	}
      }
 
//       for(ActiveSetIterator s=_active.begin(); s!=_active.end(); s++) {
// 	DefNode *n = s->first;
// 	double f3 = s->second;
// 	double x3 = n->getPoint(2);
// 	double dx3 = x3 - _Z;
// 	_penetration += dx3*dx3;
// 	if(f0) {
// 	  _energy += f3*dx3 + 0.5*_k*dx3*dx3;
// 	}
// 	if(f1) {
// 	  f3 += _k*dx3;
// 	  _FZ += f3;
// 	  n->addForce(2,f3);
// 	}
// 	_penetration = sqrt( _penetration );
//       }
    }
  private:

    double _k;
    double _Z;
    double _FZ;
    DefNodeContainer _defNodes;
    bool _up;

//     ActiveSet _active;
    std::vector<bool> _active;
    blitz::Array<double,1> _forces;
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
      _dx = 0.0;
      _x0 = x0;
      _bc = bc;
      _f = 0.0, 0.0, 0.0;
    }

    double penetration() const {return _dx;}

    double penaltyCoefficient() const {return _k;}

    void setPenaltyCoefficient( double k ) { _k = k; }

    // E = 0.5 * k * ( x3-Z )^2
    // f3 = dE/dx3 = k * ( x3-Z )
    //! Do mechanics on Body
    virtual void compute( bool f0, bool f1, bool f2 ) {
      //std::cout << "RigidBC::compute()" << std::endl;

      if(f0) _energy = 0.0;

      Vector3D dx(0.0);
      for(int i=0; i<3; i++) {
	if(_bc(i)) {
	  dx(i) = _defNode->getPoint(i) - _x0(i);
	  _dx = std::max(std::abs(dx(i)), _dx);
	}
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
    double _dx;
    Vector3D _x0;
    Vector3D _f;
    VectorBC _bc;
    DefNode * _defNode;
  };


} // end namespace voom
