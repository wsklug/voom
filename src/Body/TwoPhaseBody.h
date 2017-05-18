// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//


/*! 
  \file TwoPhaseBody.h

  \brief TwoPhaseBody is a concrete class derived from Body.  
   The body is composed of subdivision two phase elements

*/

#if !defined(__TwoPhaseBody_h__)
#define __TwoPhaseBody_h__

#include<blitz/array.h>
#include<vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
#include "Body.h"
#include "TwoPhaseElement.h"
#include "ads/halfedge.h"
//#include "DefinedTypes.h"
#include "Node.h"
#include "HDSVertexNodeTwo.h"
#include "HDSConnectivity.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace voom
{

  template < class Material_t >
  class TwoPhaseBody : public Body
  {
  public:
    
    

    // typedefs
    typedef TwoPhaseElement<Material_t> FeElement_t;
    typedef typename FeElement_t::Node_t FeNode_t;
    typedef typename std::vector< FeNode_t* > FeNodeContainer;
    typedef typename FeNodeContainer::iterator FeNodeIterator;
    typedef typename FeNodeContainer::const_iterator ConstFeNodeIterator;

    typedef std::vector< FeElement_t* > FeElementContainer;
    typedef typename FeElementContainer::iterator FeElementIterator;
    typedef typename FeElementContainer::const_iterator ConstFeElementIterator;
		
    typedef typename FeElement_t::GlobalConstraint GlobalConstraint;

    typedef tvmet::Vector<int,3> ElementConnectivity;
    typedef std::vector<ElementConnectivity> ConnectivityContainer;

    //! Default Constructor
    TwoPhaseBody() {;}

    //! Construct from stuff
    TwoPhaseBody (Material_t material,
		  ConnectivityContainer & connectivities,
		  const NodeContainer & nodes,
		  const unsigned quadOrder = 1,
		  const int numberOfBoundaries = 0,
		  const double pressure = 0.0,
		  const double tension = 0.0,
		  const double chemicalTension = 0.0,
		  const double penaltyVolume = 1.0e4,
		  const double penaltyArea= 1.0e6,
		  const double penaltyAreaOne = 1.0e6,
		  const double viscosity=0.0,
		  GlobalConstraint volumeConstraint = FeElement_t::penalty,
		  GlobalConstraint areaConstraint = FeElement_t::penalty, 
		  GlobalConstraint areaOneConstraint = FeElement_t::penalty) {
      initializeBody(material, connectivities, nodes, quadOrder, 
		     numberOfBoundaries, pressure, tension, chemicalTension, 
		     penaltyVolume, penaltyArea, penaltyAreaOne, viscosity,
		     volumeConstraint, areaConstraint, areaOneConstraint);
    }

    //! initialize
    void initializeBody(Material_t material,
			ConnectivityContainer & connectivities,
			const NodeContainer &  nodes,
			const unsigned quadOrder,
			const int numberOfBoundaries,
			const double pressure,
			const double tension,
			const double chemicalTension,
			const double penaltyVolume,
			const double penaltyArea,
			const double penaltyAreaOne,
			const double viscosity,
			GlobalConstraint volumeConstraint,
			GlobalConstraint areaConstraint,
			GlobalConstraint areaOneConstraint );
    
    //! virtual destructor
    virtual ~TwoPhaseBody() {;}
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    double volume() { return _volume; }
    double constraintVolume() { return _constraintVolume; }
    void setConstraintVolume(double V) { _constraintVolume = V; }
    
    //! get total area
    double area() { return _area; }
    double constraintArea() { return _constraintArea; }
    void setConstraintArea(double A) { _constraintArea = A; }
    
    //! get area of the first component    
    double areaOne() { return _areaOne; }
    double constraintAreaOne() { return _constraintAreaOne; }
    void setConstraintAreaOne(double A) { _constraintAreaOne = A; }

    //! get pressure
    double pressure() const {return _pressureNode->point();}

    //! get tension
    double tension() const {return _tensionNode->point();}

    //! get chemical tension
    double chemicalTension() const {return _chemicalTensionNode->point();}

    // Copy field values from nodes into arrays
        
    //! set constraint volume
    void reduceVolume(const double nu) {
      double A = _constraintArea;
      _constraintVolume = nu*sqrt(A*A*A/M_PI)/6.0;
    }

    void setCytoSpring(double new_mu, double new_kS, double new_kSpring){
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	(*e)->setCytoSpring(new_mu, new_kS, new_kSpring);
      }	      
    }

    //
    //
    //
    void printWork() const{
      std::cout<<"Energy by volume: "<< _tempWork1 << std::endl
	       <<"Energy by area  : "<< _tempWork2 << std::endl;
    }
    //! loop over elements and add up the constraint energy of each one
    double constraintEnergy() const {
      double constraintEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	constraintEnergy += (*e)->localConstraintEnergy();
      }	
      return constraintEnergy;
    }

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
	
      }	
      return totalStrainEnergy;
    }
    
    //! loop over elements and add up the work of each one
    //double work() const {
    //double work = 0.0;
    //for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
    //work += -(*e)->work()+(*e)->constraintEnergy()+(*e)->chemicalConstraintEnergy();
    //}	
    //return work;
    //}

    //! print body connectivities by using Half-edge Data Structure info.
    void printByHDS();
        
    //! output to openDX format
    void createOpenDXData(const std::string& fileName, 
			  const int nWhichForce = 2);
    void printParaview(const std::string fileName) const ;
    
    //! recreate input file
    void createInputFile(const std::string& filename);
    
    virtual void incrementLoads(const double load){}

    void resetReference() {
      for(FeNodeIterator n = _shellNodes.begin(); n!= _shellNodes.end(); n++) {
	(*n)->resetPosition();
      }

      for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	(*e)->updateRefConfiguration();
      }	
      return;
    }
    
    // access viscous regularization "energy"
    double viscousEnergy() const {return _viscousEnergy;}

    void setViscosity (const double & vis){
      _viscosity = vis;
    }

    //
    // Compute force along the generatix in the vescile (Paul's model)
    // for comparison
//     void ComputeVescileForce();
    
  private:

    // typedefs from namespace ads for halfedge data structure
    
    //! Half-edge data structure
    typedef ads::HalfedgeDS<HDSVertexNode, ads::HDSHalfedge, HDSConnectivity> HDS;
    
    //! Iterator for vertices of a HDS 
    typedef HDS::Vertex_iterator Vertex_iterator;
    
    //! Handle to a vertex
    typedef HDS::Vertex_handle Vertex_handle;
    
    //! const handle to a vertex
    typedef HDS::Vertex_const_handle Vertex_const_handle;
    
    //! Handle to a halfedge
    typedef HDS::Halfedge_handle Halfedge_handle;
    
    //! const handle to a halfedge
    typedef HDS::Halfedge_const_handle  Halfedge_const_handle;
    
    //! handle to a face
    typedef HDS::Face_handle Face_handle;
    
    //! const handle to a halfedge
    typedef HDS::Face_const_handle Face_const_handle;
    
    //! mapping from two vertices to a halfedge handle
    typedef std::map<std::pair<int, int>, Halfedge_handle> MapHalfedge;
    
    //! mapping from a vertex to a vertex handle
    typedef std::map<int, Vertex_handle> MapVertex;
    
    //! two nodal indices -> one edge
    typedef std::vector<std::pair<int, int> > IntPairsContainer;
    typedef IntPairsContainer::iterator IntPairsIterator;
    //
    //! Face handle container
    typedef std::vector < Face_handle > FaceHandleContainer;	
    typedef FaceHandleContainer::iterator FaceHandleIterator;
    //
    //! Corner Valences
    typedef tvmet::Vector<unsigned int, 3> CornerValences;
        
    //
    // data
    //
    
    // half edge data structure
    HDS _hds;

    //! Nodes
    FeNodeContainer _shellNodes;

    //! Elements
    FeElementContainer 	_shells;		


    // for volume constraint
    GlobalConstraint _volumeConstraint;
    double _volume;
    double _constraintVolume;
    double _penaltyVolume;
    MultiplierNode * _pressureNode;

    // for areaOne constraint
    GlobalConstraint _areaOneConstraint;
    double _areaOne;
    double _constraintAreaOne;
    double _penaltyAreaOne;
    MultiplierNode * _chemicalTensionNode;

    // for total area constraint
    GlobalConstraint _areaConstraint;
    double _area;
    double _constraintArea;
    double _penaltyArea;
    double _viscosity;
    double _viscousEnergy;
    MultiplierNode * _tensionNode;

    //
    //   private methods
    //
    double _tempWork1, _tempWork2;    

    void _convertConnectivityToPairs(const tvmet::Vector<int, 3>& v,
				     IntPairsContainer& ipc);
    void _reversePairs(IntPairsContainer& ipc);

    void _initializeMap(MapHalfedge& mh, 
			MapVertex& mv, 
			FaceHandleContainer& uf);

    void _createHDS(ConnectivityContainer & connectivities, 
		    const int numberOfBoundaries = 0);

    int  _newNodalIndex(const MapHalfedge::iterator itrMapHandle, 
			const tvmet::Vector<int, 3>& c);

    void _createElement(Material_t material, 
			const Face_handle fh, 
			const unsigned quadOrder = 1);

//     //! initialize nodes' neighbors
//     void _initNodeNeighbors();
    
  };  
} // namespace voom

#include "TwoPhaseBody.cc"
#include "TwoPhaseBody-HDS.cc"

#endif // __TwoPhaseBody_h__
