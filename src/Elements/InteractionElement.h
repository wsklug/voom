/*!
  \file InteractionElement.h

  \InteractionElement is derived from element. It has two LoopShell
  element objects and deals the interaction of these two LoopShell
  elements

*/
#if !defined(__InteractionElement_h__)
#define __InteractionElement_h__

#include <blitz/array.h>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include <vector>

#include "Element.h"
#include "LoopShell.h"
#include "LoopShellShape.h"
#include "Node.h"
#include "ShellGeometry.h"
#include "VoomMath.h"


namespace voom
{

  template<class Material_t>
  class InteractionElement : public Element
  {

  public:

    typedef tvmet::Vector<double, 3> Vector3D;

    typedef typename LoopShell<Material_t>::NodeContainer           NodeContainer;
    typedef typename LoopShell<Material_t>::NodeIterator            NodeIterator;
    typedef typename LoopShell<Material_t>::ConstNodeIterator       ConstNodeIterator;
    typedef typename LoopShell<Material_t>::QuadPointIterator       QuadPointIterator;
    typedef typename LoopShell<Material_t>::ConstQuadPointIterator  ConstQuadPointIterator;

    //! virtual destructor
    virtual ~InteractionElement() {;}

    //! constructor
    InteractionElement(const double k,
		       const double rZero,
		       const double epsilon,
		       LoopShell<Material_t>* e1,
		       LoopShell<Material_t>* e2
		       )
      {
	_e1=e1;
	_e2=e2;
	_k=k;
	_rZero=rZero;
	_epsilon=epsilon;
	_rc=_rZero + sqrt(_epsilon/_k);

	_nNodes1 = _e1->_nodes.size();
	_nNodes2 = _e2->_nodes.size();

	_internalForce1.resize( _nNodes1 );
	_internalForce2.resize( _nNodes2 );

	//store nodes in baseNodes for consistency check of Element
	for(ConstNodeIterator n=(_e1->_nodes).begin(); n!=(_e1->_nodes).end(); n++)
	  _baseNodes.push_back(*n);
	for(ConstNodeIterator n=_e2->_nodes.begin(); n!=_e2->_nodes.end(); n++)
	  _baseNodes.push_back(*n);


      }

    double adhesionEnergy() const { return _adhesionEnergy; }

    virtual  void compute(bool f0, bool f1, bool f2);

    double distance() const {
      DeformationNode<3>::Point x1avg(0.0), x2avg(0.0);
      for (int a = 0; a < 3; a++){
	x1avg+= _e1->_nodes[a]->point();
	x2avg+= _e2->_nodes[a]->point();
      }

      x1avg/=3.0;
      x2avg/=3.0;
      double distance = tvmet::norm2(x1avg - x2avg);
      return distance;
    }

    //if the two elements are not connected
    bool connectionFlag() const {    
      for (int a = 0; a < _e1->_nodes.size(); a++){
	int x1a = _e1->_nodes[a]->id();
	
	for(int b = 0; b < _e2->_nodes.size(); b++){
	  int x2b = _e2->_nodes[b]->id();
	  if(x1a == x2b)
	    return true;
	}
      }
      return false;
    }

  private:

    double _k;
    double _rZero;
    double _rc;
    double _epsilon;
    int    _nNodes1;
    int    _nNodes2;

    double _adhesionEnergy;

    blitz::Array< Vector3D, 1> _internalForce1; //forces on e1
    blitz::Array< Vector3D, 1> _internalForce2; //forces on e2

    LoopShell<Material_t>* _e1;
    LoopShell<Material_t>* _e2;

  };//end of class

} // namespace voom

#include "InteractionElement.cc"

#endif // __InteractionElement_h__



