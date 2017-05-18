// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    HoHai Van and William S. Klug
//                University of California Los Angeles
//                    (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__CardiacPotential_h__)
#define __CardiacPotential_h__

#include <vector>
#include "Node.h"
#include "NodeBase.h"
#include "Element.h"


namespace voom
{
  //! A position node with scalar point for voltage, also storing damping coefficients.
  template< int dim_n >
  class VoltageNode 
    : public PositionNode< double, dim_n > {
    
  public:
    typedef PositionNode< double, dim_n > Base;
    typedef double Point;
    typedef typename tvmet::Vector<double,dim_n> PositionVector;

//     using Node<Point>::_id;
//     using Node<Point>::_point;
//     using Node<Point>::_force;

    using Base::_id;
    using Base::_point;
    using Base::_force;
    using Base::_X;



    //! default constructor
    VoltageNode() {}

    //! construct from position and point
    VoltageNode(int id, const NodeBase::DofIndexMap & index, 
		const PositionVector & X, const Point p) : Base(id,index)
    { _X = X; _point = p; _damping = 0.0;}

    // it seems like this should be inherited from Node, but GCC objects...?
    virtual void setPoint( const Point & p ) { _point = p; }


    double getPoint(int i) const { assert(i<1); return _point; }

    void setPoint(int i, double x) {assert(i<1); _point = x; }
   
    void addPoint(int i, double dx) {assert(i==0); _point += dx;}

    double getForce(int i) const { assert(i<1); return _force; }

    void setForce(int i, double f) {assert(i<1); _force = f; }

    void addForce(int i, double df) {assert(i==0); _force += df; }

    virtual void setForce( const Point & f ) { _force = f; }

    double getDamping(int i) const { assert(i<1); return _damping; }

    void setDamping(int i, double d) { assert(i<1); _damping = d; }

    void addDamping(int i, double d) {assert (i==0);_damping +=d; }

    int dof() const {return 1;}

  protected:

    double _damping;

  };


  //! CardiacPotential is an element for cardiac electrical propagation.
  /*! This class implements an Element for cardiac electrical
    propagation by adding vertex nodes which store voltage DOF as well
    as nodes at the quadrature points which store internal DOF.  Shape
    function and material objects are stored inside quadrature point
    structs.  This allows for convenient iteration through a container
    of quad points as is typically done for computation of element
    integrals.  The Quadrature_t object is to be used in initializing
    the weights and positions of the quad points.
   */

  template< class Quadrature_t,
	    class Shape_t,
            const int dim_t >
  class CardiacPotential 
  {
    
  public:

    /*! VoltageNodes are placed at the element vertices for
     *  isoparametric interpolation of the (scalar) voltage field.
     *  Nodal voltages are stored in the data member _point of the
     *  position node class.
     */
// This is for 2D    
    typedef VoltageNode< 2 > Node_t; 

// This is for 3D
//  typedef VoltageNode<3> Node_t;

    typedef  std::vector<Node_t*> VoltageNodeContainer;
    typedef  VoltageNodeContainer::iterator VoltageNodeIterator;
    typedef  VoltageNodeContainer::const_iterator ConstVoltageNodeIterator;

    //! Nodes for internal DOF (relating to ionic concentrations, etc.)
    //! 7 DOF for each gating variable
    typedef VectorNode< 9 > InternalNode;
    typedef std::vector<InternalNode*> InternalNodeContainer;
    typedef InternalNodeContainer::iterator InternalNodeIterator; 
    typedef InternalNodeContainer::const_iterator ConstInternalNodeIterator;



    typedef blitz::Array<double, 2> StiffnessMatrix;

    //! default constructor
    CardiacPotential( const Quadrature_t & quad,
		      const VoltageNodeContainer & nodes
		     );

    //! virtual destructor
    virtual ~CardiacPotential() {;}
    
    //! Do electrophysiology on element; compute current, voltage, etc.
    virtual void compute_v(double istim);
    virtual void compute_gate(double dt);

    virtual void computeLumpedCapacitance();


    //! Access the container of nodes
    const VoltageNodeContainer& nodes() const { return _vNodes; }

    //! Access the container of nodes
    const VoltageNodeContainer& voltageNodes() const { return _vNodes; }

    //! Access the container of nodes
    InternalNodeContainer& internalNodes() { return _iNodes; }

    //! Structure to hold together everything needed at a quad point
    /*! Stuff needed at each quad point includes the quad weight,
      shape functions and their derivatives evaluated at the point,
      and an internal node which stores the internal DOF at the quad
      point.
     */
    struct QuadPointStruct {
      double weight;
      typename Shape_t::FunctionContainer shapeFunctions;
      typename Shape_t::DerivativeContainer shapeDerivatives;
      InternalNode* internalNode;
      QuadPointStruct(double w, const Shape_t & s, 
		      const VoltageNodeContainer & nds);
      
       };

    typedef typename std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
   		    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }
    QuadPointContainer & quadraturePoints() { return _quadPoints; }

    
  protected:
    VoltageNodeContainer _vNodes;  
    InternalNodeContainer _iNodes;  
    QuadPointContainer _quadPoints;
    StiffnessMatrix _stiffness;


  };
	
} // namespace voom

#endif // __CardiacPotential_h__
