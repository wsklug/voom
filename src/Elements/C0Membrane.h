// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file C0Membrane.h

  \brief C0Membrane is a concrete class derived from element, implementing
  a C0 finite element for finite deformation membrane structures.

*/

#if !defined(__C0Membrane_h__)
#define __C0Membrane_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "StandardElement.h"
#include "Node.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Concrete class for a C0 finite deformation membrane finite element.
  */
  template<class Quadrature_t, class Material_t, class Shape_t>
  class C0Membrane 
    : public StandardElement< DeformationNode<3>, Quadrature_t,
			      Material_t, Shape_t > 
  {
    
  public:

    typedef DeformationNode<3> Node_t;

    // WSK & MMG don't think this typedef should be here!!!
//     typedef TriangleQuadrature Quadrature_t;

    typedef StandardElement<Node_t, Quadrature_t, Material_t, Shape_t> Base;

    typedef typename Base::NodeContainer 		NodeContainer;
    typedef typename Base::NodeIterator 		NodeIterator;
    typedef typename Base::ConstNodeIterator 		ConstNodeIterator;
    typedef typename Base::QuadPointIterator 		QuadPointIterator;
    typedef typename Base::ConstQuadPointIterator 	ConstQuadPointIterator;

    using Base::_nodes;
    using Base::_quadPoints;
    using Base::_baseNodes;
    using Base::_energy;

    //! virtual destructor
    virtual ~C0Membrane() {;}

    //! Constructor
    C0Membrane( const Quadrature_t & quad,
		const Material_t & mat,
		const NodeContainer & nodes,
		MultiplierNode * pressureNode,
		MultiplierNode * tensionNode,
		GlobalConstraint volumeConstraint = noConstraint,
		GlobalConstraint areaConstraint = noConstraint
	       )
    {
      //! initialize NodeContainer
      unsigned nNodes = nodes.size();
      _nodes = nodes;

      for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	_baseNodes.push_back(*n);
      
      _pressureNode = pressureNode;
      _tensionNode = tensionNode;

      _areaConstraint = areaConstraint;
      _volumeConstraint = volumeConstraint; 

      _volume = 0.0;
      _area = 0.0;
	  
      //! initialize materials and shape functions
      _quadPoints.clear();
      for(typename Quadrature_t::ConstPointIterator p=quad.begin(); 
	  p!=quad.end(); p++) {
	Shape_t shp( p->coords );
	_quadPoints.push_back( typename Base::QuadPointStruct(p->weight, mat, shp) ); 
	if( shp.functions().size() != nNodes ) {
	  std::cout << "Number of nodes: " << nNodes << std::endl
		    << "Number of functions: " << shp.functions().size()
		    << std::endl
		    << "These should be equal." << std::endl;
	  exit(0);
	}
      }

      //! allocate memory for mechanics variables
      _internalForce.resize( nNodes );
      _pressureForce.resize( nNodes );
      _tensionForce.resize( nNodes  );

      updateRefConfiguration(); 

      //! initialize values = 0 by default
      //compute(false, false, false); 
      compute(true, true, false); 
    }


  public:

//     void setMaterialY(double Y){
//       for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
// 	Material_t& material = p->material;
// 	material.setYoungsModulus(Y);
// 	//material.setCytoSpring(0.0,0.0,kSpring);
//       }
//     }
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }

    //! Accessor for element work done by pressure
    double work() const { return _work; }
		
    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }

    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration();

    //
    //   data
    //
  private:
    //! forces on the element
    blitz::Array< Vector3D, 1> _internalForce;
    blitz::Array< Vector3D, 1> _pressureForce;
    blitz::Array< Vector3D, 1> _tensionForce;
		
    //! uniform pressure applied on the element
    MultiplierNode * _pressureNode;

    MultiplierNode * _tensionNode;

    //! enum for area constraint choice
    GlobalConstraint _areaConstraint;
    //! enum for volume constraint choice
    GlobalConstraint _volumeConstraint;

    //! strain energy
    double _strainEnergy;

    //! work done by pressure
    double _work;

    //! volume corresponding to this element
    double _volume;

    //! area enclosed by the current element
    double _area;

  };  // end of class
} // namespace voom

#include "C0Membrane.icc"

#endif // __C0Membrane_h__
