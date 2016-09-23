// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------
// 

/*! 
  \file C1SplineAxiShell.h

  \brief C1SplineAxiShell is a concrete class derived from element, implementing
  a thin shell subdivision surface finite element.

*/

#if !defined(__C1SplineAxiShell_h__)
#define __C1SplineAxiShell_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "StandardElement.h"
#include "LineQuadrature.h"
#include "Spline.h"
#include "Node.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Concrete class for an axisymmetric shell element with C1
    interpolation for the position field. Node 1 and 2 form the
    element, with the supporting region extending to the neighbouring
    node 0 and 3
  */
  template<class Material_t>
  class C1SplineAxiShell 
    : public StandardElement< DeformationNode<2>, LineQuadrature,
			      Material_t, Spline > 
  {
    
  public:

    typedef DeformationNode<2> Node_t;
    typedef Node_t::Point Point_t;

    typedef LineQuadrature Quadrature_t;
    typedef Spline Shape_t;

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
    virtual ~C1SplineAxiShell() {;}

    //! Constructor
    C1SplineAxiShell( const Quadrature_t & quad,
	       const Material_t & mat,
	       const NodeContainer & nodes,
	       MultiplierNode * pressureNode = 0,
	       MultiplierNode * tensionNode = 0,
               MultiplierNode * totalCurvatureNode = 0,
	       GlobalConstraint volumeConstraint = noConstraint,
	       GlobalConstraint areaConstraint = noConstraint,
               GlobalConstraint totalCurvatureConstraint = noConstraint
	       )

  {
      //! initialize NodeContainer
      unsigned nNodes = nodes.size();
      assert( nNodes == 4 );

      _nodes = nodes;

      for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
	_baseNodes.push_back(*n);
      
      _pressureNode = pressureNode;
      _tensionNode = tensionNode;
      _totalCurvatureNode = totalCurvatureNode;
      _areaConstraint = areaConstraint;
      _volumeConstraint = volumeConstraint;   
      _totalCurvatureConstraint = totalCurvatureConstraint;   

      _volume = 0.0;
      _area = 0.0;
      _totalCurvature = 0.0;
	  
      //! initialize materials and shape functions
      _quadPoints.clear();
      for(typename Quadrature_t::ConstPointIterator p=quad.begin(); 
	  p!=quad.end(); p++) {
	Shape_t shp( p->coords[0] );
	_quadPoints.push_back( typename Base::QuadPointStruct(p->weight, mat, shp) ); 
      }

      //! allocate memory for mechanics variables
      _internalForce.resize( nNodes );
      _pressureForce.resize( nNodes );
      _tensionForce.resize( nNodes );

      updateRefConfiguration(); 
      compute(false,false,false);
      _prescribedVolume = _volume;
      _prescribedArea = _area;
      _prescribedTotalCurvature = _totalCurvature;

      //! initialize values = 0 by default
      compute(true, true, false); 
    }


  public:
   
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    virtual void compute(bool f0, bool f1, bool f2);
    
    double strainEnergy() const { return _strainEnergy; }

    double bendingEnergy() const { return _bendingEnergy; }

    double stretchingEnergy() const { return _stretchingEnergy; }

    //! Access the container of nodes
    const NodeContainer & nodes() { return _nodes; }
    
    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }
    //! access total curvature
    double totalCurvature() const {return _totalCurvature;}

    //! access prescribed volume
    double prescribedVolume() const { return _prescribedVolume; }
    //! access prescribed area
    double prescribedArea() const { return _prescribedArea; }
    //! access prescribed total curvature
    double prescribedTotalCurvature() const {return _prescribedTotalCurvature;}

    //! modify prescribed volume
    double setPrescribedVolume(double V) {_prescribedVolume = V;}
    //! modify prescribed volume
    double setPrescribedArea(double A) {_prescribedArea = A;}

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

    //! uniform tension work-conjugate to global area
    MultiplierNode * _tensionNode;

    MultiplierNode * _totalCurvatureNode;

    //! enum for area constraint choice
    GlobalConstraint _areaConstraint;
    //! enum for volume constraint choice
    GlobalConstraint _volumeConstraint;

    GlobalConstraint _totalCurvatureConstraint;

    //! strain energy
    double _strainEnergy;

    //! bending energy
    double _bendingEnergy;

    //! stretching energy
    double _stretchingEnergy;

    //! work done by pressure
    double _work;

    //! volume corresponding to this element
    double _volume;

    //! prescribed volume corresponding to this element
    double _prescribedVolume;

    //! area enclosed by the current element
    double _area;

    //! prescribed area enclosed by the current element
    double _prescribedArea;

    //! integral of mean curvature
    double _totalCurvature;

    //! prescribed integral of mean curvature
    double _prescribedTotalCurvature;
	
  };  // end of class
} // namespace voom

#include "C1SplineAxiShell.icc"

#endif // __C1SplineAxiShell_h__
