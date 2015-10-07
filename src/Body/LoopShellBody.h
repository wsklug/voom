// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file LoopShellBody.h

  \brief LoopShellBody is a concrete class derived from Body, implementing
  the concept of a thin shell body.  The body is composed of
  subdivision shell elements

*/

#if !defined(__LoopShellBody_h__)
#define __LoopShellBody_h__

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
#include "LoopShell.h"
#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "HalfEdgeMesh.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{


  /*!  
	Concrete class for a body composed of subdivision shell elements
  */
  template < class Material_t >
  class LoopShellBody : public Body
  {
  public:
    
    

    // typedefs
    typedef LoopShell<Material_t> FeElement_t;
    typedef typename FeElement_t::Node_t FeNode_t;
    typedef typename std::vector< FeNode_t* > FeNodeContainer;
    typedef typename FeNodeContainer::iterator FeNodeIterator;
    typedef typename FeNodeContainer::const_iterator ConstFeNodeIterator;

    typedef std::vector< FeElement_t* > FeElementContainer;
    typedef typename FeElementContainer::iterator FeElementIterator;
    typedef typename FeElementContainer::const_iterator ConstFeElementIterator;
		
    typedef std::vector< Constraint* > ConstraintContainer;
    typedef typename ConstraintContainer::iterator ConstraintIterator;
    typedef typename ConstraintContainer::const_iterator ConstConstraintIterator;

    typedef tvmet::Vector<int,3> ElementConnectivity;
    typedef std::vector<ElementConnectivity> ConnectivityContainer;

    //! Default Constructor
    LoopShellBody() {;}

    //! Construct from stuff
    LoopShellBody(Material_t material,
		  ConnectivityContainer & connectivities,
		  const NodeContainer & nodes,
		  const unsigned quadOrder = 1,
		  const double pressure = 0.0,
		  const double tension = 0.0,
		  const double totalCurvatureForce = 0.0,
		  const double penaltyVolume = 1.0e4,
		  const double penaltyArea= 1.0e6,
		  const double penaltyTotalCurvature= 1.0e4,
		  GlobalConstraint volumeConstraint = noConstraint,
		  GlobalConstraint areaConstraint = noConstraint,
		  GlobalConstraint totalCurvatureConstraint = noConstraint ) {

      initializeBody(material, connectivities, nodes, quadOrder, 
		     pressure, tension, totalCurvatureForce,
		     penaltyVolume, penaltyArea, penaltyTotalCurvature, 
		     volumeConstraint, areaConstraint, totalCurvatureConstraint);
    }

    //! initialize
    void initializeBody(Material_t material,
			ConnectivityContainer connectivities,
			const NodeContainer &  nodes,
			const unsigned quadOrder,
			const double pressure,
			const double tension,
			const double totalCurvatureForce,
			const double penaltyVolume,
			const double penaltyArea,
			const double penaltyTotalCurvature,
			GlobalConstraint volumeConstraint,
			GlobalConstraint areaConstraint,
			GlobalConstraint totalCurvatureConstraint);

    
    //! virtual destructor
    virtual ~LoopShellBody() {
      delete _pressureNode;
      delete _tensionNode;
      delete _totalCurvatureNode;
      for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	delete (*e);
      }
    }

    FeElementContainer & shells() {return _shells;};
    FeNodeContainer & shellsNodes() {return _shellNodes;};
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );

    //! calculate the curvatures at all elements
    void cal_curv(std::vector<double> &curv);
    
    double volume() const{ return _volume; }
    double prescribedVolume() const { return _prescribedVolume; }
    void setPrescribedVolume(double V) { _prescribedVolume = V; }
    
    //! get total area
    double area() const { return _area; }
    double prescribedArea() const { return _prescribedArea; }
    void setPrescribedArea(double A) { _prescribedArea = A; }

    double totalCurvature() const { return _totalCurvature;}
    double prescribedTotalCurvature() const { return _prescribedTotalCurvature;}
    
    void SetRefConfiguration(double edge);

    //! get pressure
    double pressure() const {return _pressureNode->point();}

    //! get tension
    double tension() const {return _tensionNode->point();}

    // Copy field values from nodes into arrays
        
    //! set prescribed volume
    void reduceVolume(const double nu) {
      double A = _prescribedArea;
      _prescribedVolume = nu*sqrt(A*A*A/M_PI)/6.0;
//       double volumeFactor = nu*sqrt(A*A*A/M_PI)/6.0/_prescribedVolume;
//       _prescribedVolume = 0.0;
//       for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
// 	(*e)->setPrescribedVolume( volumeFactor*(*e)->prescribedVolume() );
// 	_prescribedVolume += (*e)->prescribedVolume();
//       }

    }

    //! set prescribed total curvature
    void setCTC (const double CTC) {
      _prescribedTotalCurvature = CTC;
    }

    //! loop over elements and add up the stretching energy of each one
    double strechingEnergy() const {
      double stretchingEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	stretchingEnergy += (*e)->strechingEnergy();
      }	
      return strechingEnergy;
    }

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
	
      }	
      return totalStrainEnergy;
    }


    //! print body connectivities by using Half-edge Data Structure info.
    void printByHDS();
        
    //! output to openDX format
    void createOpenDXData(const std::string& fileName, 
			  const int nWhichForce = 2);

    void printParaview(const std::string fileName) const;

    void printObj(const std::string name) const;
    
    //    virtual void incrementLoads(const double load){}

    void resetReference() {
      for(FeNodeIterator n = _shellNodes.begin(); n!= _shellNodes.end(); n++) {
	(*n)->resetPosition();
      }

//       for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
// 	(*e)->updateRefConfiguration();
//       }	
//      cytoskeleton energy not updated

      //compute(true,true,true);//add this because we have nodal force add up twice (updateRefConfiguration() has compute)
                              //zero out force first. 
      return;
    }

    virtual void pushBack( Element* e ) { 
      _elements.push_back(e); 
      _active.push_back(true);}

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    double fixedPressure() const { return _fixedPressure; }

    double fixedTension() const { return _fixedTension; }

    void updateFixedPressure() {
      _fixedPressure = _pressureNode->point();
    }

    void updateFixedTension() {
      _fixedTension = _tensionNode->point();
    }


    void updateFixedForce(double P, double T, double TC){
      _fixedPressure = P;
      _fixedTension  = T;
      _fixedTotalCurvatureForce = TC;
    }

    void updatePenaltyVolumeArea(double pV, double pA, double pTC){
      _penaltyVolume = pV;
      _penaltyArea = pA;
      _penaltyTotalCurvature = pTC;
    }

    //! Query an element's activity status
    bool active(int e) { return _active[e]; }

    //! Mark an element as active so it will be computed
    void activate(int e) { _active[e] = true; }

    //! Mark an element as inactive so it will not be computed
    void deactivate(int e) { _active[e] = false; }



    // ----------------------------------------- //
    // New functions //

    // Compute average edge length
    double AverageEdgeLength();

    // Compute element neighbors
    std::vector<std::vector<uint > > ComputeElementNeighBors(std::vector<tvmet::Vector<int,3> > ConnTable);

    // Remesh elements with bad aspect ratio
    uint Remesh(double ARtol, Material_t material, uint quadOrder);

    // Create Shell FE
    void CreateLoopFE(ConnectivityContainer & connectivities, Material_t material, uint quadOrder, bool remeshing);

    // Calculate and return maximum principal strain (using right
    // Cauchy-Green strain) for all active elements
    std::vector<double> calcMaxPrincipalStrains() const;

    //Get Maximum Principal Strains
    //std::vector<double> getMaxPrincipalStrains() {return _maxPrincipalStrain;};

    void setAreaConstraint(GlobalConstraint AreaConstr) {_areaConstraint = AreaConstr;};
    void setVolumeConstraint(GlobalConstraint VolConstr) {_volumeConstraint = VolConstr;};
    // ----------------------------------------- //



  private:

    //! Store corner valences in a vector for building irregular elements
    typedef tvmet::Vector<unsigned int, 3> CornerValences;
        
    //! Nodes
    FeNodeContainer _shellNodes;

    //! Elements
    FeElementContainer 	_shells;		

    //! vector of boolean flags to activate/deactivate elements
    std::vector<bool> _active;

    //! vector of largest eigen value of right Cauchy Green strain for
    //! each element
    //std::vector<double> _maxPrincipalStrain;

    // for volume constraint
    GlobalConstraint _volumeConstraint;
    double _volume;
    double _prescribedVolume;
    double _penaltyVolume;
    MultiplierNode * _pressureNode;

    // for total area constraint
    GlobalConstraint _areaConstraint;
    double _area;
    double _prescribedArea;
    double _penaltyArea;
    MultiplierNode * _tensionNode;

    // for total curvature constraint
    GlobalConstraint _totalCurvatureConstraint;
    double _totalCurvature;
    double _prescribedTotalCurvature;
    double _penaltyTotalCurvature;
    MultiplierNode * _totalCurvatureNode;

    ConstraintContainer _constraints;

    double _fixedPressure, _fixedTension, _fixedTotalCurvatureForce;

    //std::vector<int> _dimpleElementsList;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
    // Compute the diagonal elements of the stiffness matrix by
    // numerical differentiation
    void _computeStiffness();
    
  };  
} // namespace voom

#include "LoopShellBody.cc"

#endif // __LoopShellBody_h__
