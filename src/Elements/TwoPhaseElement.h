// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------


#if !defined(__TwoPhaseElement_h__)
#define __TwoPhaseElement_h__


#include <blitz/array.h>
 
 
#include <vector>
#include <cstdio>
#include <ctime>
#include <iostream>

#include "TriangleQuadrature.h"
#include "LoopShellShape.h"
#include "Node.h"
#include "Element.h"
#include "ShellGeometry.h"
#include "VoomMath.h"


namespace voom
{

  /*TwoPhaseElement defines a new element from scratch
   */

  template <class Material_t>
  class TwoPhaseElement : public Element
  {
    
  public:

    typedef tvmet::Vector<double, 4> Vector4D;

    typedef XCNode<3> Node_t;        
    typedef TriangleQuadrature Quadrature_t;
    typedef LoopShellShape Shape_t;

    enum GlobalConstraint {
      none = 0,
      multiplier = 1,
      penalty = 2
    };
    
    //! default constructor
    TwoPhaseElement() {;}

    //! destructor
    ~TwoPhaseElement() {;}

    typedef typename std::vector<Node_t*> NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;   

    using Element::_baseNodes;
    using Element::_energy;

    struct QuadPointStruct 
    {
      double	         weight;
      Material_t         material;
      LoopShellShape 	 shape;

      QuadPointStruct(double w, Material_t m, Shape_t s) 
	: weight(w), material(m), shape(s) {}
    };

    typedef std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
   		    
    TwoPhaseElement( const Quadrature_t & quad,
		     const Material_t & mat,
		     const NodeContainer & nodes,
		     const Shape_t::CornerValences& V,
		     MultiplierNode * pressureNode,
		     MultiplierNode * tensionNode,
		     MultiplierNode * chemicalTensionNode,
		     GlobalConstraint volumeConstraint = penalty,
		     GlobalConstraint areaConstraint = penalty,
		     GlobalConstraint areaOneConstraint = penalty
		     )
      {
	//! initialize NodeContainer
	unsigned nNodes = nodes.size();
	assert(nNodes == V(0)+V(1)+V(2) - 6);
	_nodes = nodes;

	for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	  _baseNodes.push_back(*n);

	_pressureNode = pressureNode;
	_tensionNode = tensionNode;
	_chemicalTensionNode = chemicalTensionNode;
	_areaConstraint = areaConstraint;
	_volumeConstraint = volumeConstraint;      
	_areaOneConstraint = areaOneConstraint;

	_volume = 0.0;
	_area = 0.0;
	_areaOne = 0.0;
	  
	//! initialize materials and shape functions
	_quadPoints.clear();
	for(typename Quadrature_t::ConstPointIterator p=quad.begin(); 
	    p!=quad.end(); p++) {
	  Shape_t shp( nNodes, V, p->coords );
	  _quadPoints.push_back( QuadPointStruct(p->weight, mat, shp) ); 
	}

	//! allocate memory for mechanics variables
	_internalForce.resize( nNodes );
	_chemicalTensionForce.resize( nNodes );
	_pressureForce.resize( nNodes );
	_tensionForce.resize( nNodes );

	updateRefConfiguration(); //updateRefConfiguration() has compute(true, true, ture) 
	_constraintVolume = _volume;
	_constraintArea = _area;
	
      }
  
  public:

    double energy() const { return _energy; } 
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    void compute(bool f0, bool f1, bool f2);

    double strainEnergy() const { return _strainEnergy; }

    double localConstraintEnergy() const { return _localConstraintEnergy; }		
    
    double constraintEnergy() const { return _constraintEnergy; }

    double chemicalConstraintEnergy() const { return _chemicalConstraintEnergy; }

    double work() const { return _work; }

    double volume() const { return _volume; }

    double area() const { return _area; }

    double areaOne() const { return _areaOne; }

    double constraintVolume() const { return _constraintVolume; }

    double constraintArea() const { return _constraintArea; }

    double constraintAreaOne() const { return _constraintAreaOne; }

    double setConstraintVolume(double V) {_constraintVolume = V;}

    double setConstraintArea(double A) {_constraintArea = A;}

    double setConstraintAreaOne(double A1) {_constraintArea = A1;}
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    void setCytoSpring(double new_mu, double new_kS, double new_kSpring){
      for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
	Material_t& material = p->material;
	material.setCytoSpring(new_mu, new_kS, new_kSpring);
      }
    }

    //! Access the container of nodes
    const NodeContainer& nodes() const { return _nodes; }

    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }  

  protected:

    NodeContainer _nodes;  
    QuadPointContainer _quadPoints;

    //! forces on the element
    blitz::Array< Vector4D, 1> _internalForce;
    blitz::Array< Vector4D, 1> _chemicalTensionForce;
    blitz::Array< Vector3D, 1> _pressureForce; 
    blitz::Array< Vector3D, 1> _tensionForce;  
		
    //! uniform pressure applied on the element
    MultiplierNode * _pressureNode;

    //! uniform tension work-conjugate to global area
    MultiplierNode * _tensionNode;

    MultiplierNode * _chemicalTensionNode;

    //! enum for area constraint choice
    GlobalConstraint _areaConstraint;

    //! enum for volume constraint choice
    GlobalConstraint _volumeConstraint;

    GlobalConstraint _areaOneConstraint;

    //total energy
    //double _energy;  

    //! strain energy
    double _strainEnergy;

    //! local constraint energy (by Cytoskeleton)
    double _localConstraintEnergy;

    //! energy caused by constraint on Area 
    double _constraintEnergy;

    double _chemicalConstraintEnergy;

    double _springEnergy;

    //! work done by pressure (constraint on volume)
    double _work;

    //! volume corresponding to this element
    double _volume;
    double _constraintVolume; //related to reference scheme

    //! area enclosed by the current element
    double _area;
    double _constraintArea; //related to reference scheme

    double _areaOne;
    double _constraintAreaOne; //related to reference scheme

    double _kSpring; //spring constant of each side of triagular element

  };
	
} // namespace voom

#include "TwoPhaseElement.cc"

#endif // __TwoPhaseElement_h__
