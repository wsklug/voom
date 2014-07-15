// -*- C++ -*-
//----------------------------------------------------------------------
//
//                 William S. Klug, Feng Feng, Lin Ma
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file LoopShell.h

  \brief LoopShell is a concrete class derived from element, implementing
  a thin shell subdivision surface finite element.

*/

#if !defined(__LoopShell_h__)
#define __LoopShell_h__

#include <blitz/array.h>
 
 
#include <vector>
#include <cstdio>
#include <ctime>

#include "StandardElement.h"
#include "TriangleQuadrature.h"
#include "LoopShellShape.h"
#include "Node.h"
#include "ShellGeometry.h"
#include "VoomMath.h"

namespace voom
{

  /*!  Concrete class for a mixed inextensible subdivision shell element
    with interpolation for the position field and tension field.
  */
  template<class Material_t>
  class LoopShell 
    : public StandardElement< DeformationNode<3>, TriangleQuadrature,
			      Material_t, LoopShellShape > 
  {
    
  public:

    typedef DeformationNode<3> Node_t;
    typedef TriangleQuadrature Quadrature_t;
    typedef LoopShellShape Shape_t;

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
    virtual ~LoopShell() {;}

    //! Constructor
    LoopShell( const Quadrature_t & quad,
	       const Material_t & mat,
	       const NodeContainer & nodes,
	       const LoopShellShape::CornerValences& V,
	       MultiplierNode * pressureNode,
	       MultiplierNode * tensionNode,
               MultiplierNode * totalCurvatureNode,
	       GlobalConstraint volumeConstraint = noConstraint,
	       GlobalConstraint areaConstraint = noConstraint,
               GlobalConstraint totalCurvatureConstraint = noConstraint
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
	Shape_t shp( nNodes, V, p->coords );
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

    //! Accessor for element work done by pressure
    double work() const { return _work; }
		    
    //! access volume
    double volume() const { return _volume; }
    //! access area
    double area() const { return _area; }
    //! access prescribed volume
    double prescribedVolume() const { return _prescribedVolume; }
    //! access prescribed area
    double prescribedArea() const { return _prescribedArea; }

    double totalCurvature() const {return _totalCurvature;}

    void setPrescribedVolume(double V) {_prescribedVolume = V;}
    void setPrescribedArea(double A) {_prescribedArea = A;}


    //! compute positions
    Vector3D computePosition(const double s1, const double s2);

    //! check new positions
    void checkPositions();
		
    //! recompute the reference geometry and send it to materials
    void updateRefConfiguration(double edgeLength=-1.0);

//     void setCytoSpring(double new_mu, double new_kS, double new_kSpring){
//       for(QuadPointIterator p=_quadPoints.begin(); p!=_quadPoints.end(); p++){
// 	Material_t& material = p->material;
// 	material.setCytoSpring(new_mu, new_kS, new_kSpring);
//       }
//     }

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

    //! volume corresponding to this element
    double _prescribedVolume;

    //! area enclosed by the current element
    double _area;

    //! area enclosed by the current element
    double _prescribedArea;

    double _totalCurvature;

    double _prescribedTotalCurvature;

    double _kSpring;

    double _springEnergy;
	
  };  // end of class
} // namespace voom

#include "LoopShell.cc"

#endif // __LoopShell_h__
