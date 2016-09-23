// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file C1AxiShellBody.h

  \brief C1AxiShellBody is a concrete class derived from Body, implementing
  the concept of an axisymmetric thin shell body composed of
  axisymmetric Hermitian shell elements

*/

#if !defined(__C1AxiShellBody_h__)
#define __C1AxiShellBody_h__

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
#include "C1AxiShell.h"
#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{


  /*!  
	Concrete class for a body composed of subdivision shell elements
  */
  template < class Material_t >
  class C1AxiShellBody : public Body
  {
  public:
    
    

    // typedefs
    typedef C1AxiShell<Material_t> FeElement_t;
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
		
    typedef tvmet::Vector<int,4> ElementConnectivity;
    typedef std::vector<ElementConnectivity> ConnectivityContainer;

    //! Default Constructor
    C1AxiShellBody() {;}

    //! Construct from stuff
    C1AxiShellBody(Material_t material,
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
			ConnectivityContainer & connectivities,
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
    virtual ~C1AxiShellBody() {
      delete _pressureNode;
      delete _tensionNode;
      for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	delete (*e);
      }
    }
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    double volume() const{ return _volume; }
    double prescribedVolume() const { return _prescribedVolume; }
    void setPrescribedVolume(double V) { _prescribedVolume = V; }
    
    //! get total area
    double area() const { return _area; }
    double prescribedArea() const { return _prescribedArea; }
    void setPrescribedArea(double A) { _prescribedArea = A; }

    double totalCurvature() const { return _totalCurvature;}
    double prescribedTotalCurvature() const { return _prescribedTotalCurvature;}
    //! get pressure
    double pressure() const {return _pressureNode->point();}

    //! get pressure
    double fixedPressure() const {return _fixedPressure;}

    //! get tension
    double tension() const {return _tensionNode->point();}

    //! get tension
    double fixedTension() const {return _fixedTension;}

    //! set tension
    void setFixedTension(double t)  {_fixedTension=t;}

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

    //! loop over elements and add up the constraint energy of each one
    double stretchingEnergy() const {
      double stretchingEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	stretchingEnergy += (*e)->stretchingEnergy();
      }	
      return stretchingEnergy;
    }

    //! loop over elements and add up the strain energy of each one
    double bendingEnergy() const {
      double bendingEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	bendingEnergy += (*e)->bendingEnergy();
	
      }	
      return bendingEnergy;
    }
    
    //! loop over elements and add up the work of each one
    //double work() const {
    //double work = 0.0;
    //for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
    //work += -(*e)->work()+(*e)->constraintEnergy();
    //}	
    //return work;
    //}


    void printParaview(const std::string fileName) const ;

    virtual void incrementLoads(const double load){}

    void resetReference() {
      for(FeNodeIterator n = _shellNodes.begin(); n!= _shellNodes.end(); n++) {
	(*n)->resetPosition();
      }
      for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	(*e)->updateRefConfiguration();
      }	

      compute(true,true,true);

      return;
    }

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    void updatePenaltyVolume(double pV){
      _penaltyVolume = pV;
    }

    void updatePenaltyArea(double pA){
      _penaltyArea = pA;
    }

    void updatePenaltyCurvature(double pTC){
      _penaltyTotalCurvature = pTC;
    }

    void updateFixedPressure() { _fixedPressure = _pressureNode->point(); }

    void updateFixedTension() { _fixedTension = _tensionNode->point(); }

    void updateFixedTotalCurvatureForce() { _fixedTotalCurvatureForce = _totalCurvatureNode->point(); }

    //
    // Compute force along the generatix in the vescile (Paul's model)
    // for comparison
//     void ComputeVescileForce();
    
  private:

    //! Nodes
    FeNodeContainer _shellNodes;

    //! Elements
    FeElementContainer 	_shells;		


    // for volume constraint
    GlobalConstraint _volumeConstraint;
    double _volume;
    double _prescribedVolume;
    double _penaltyVolume;
    MultiplierNode * _pressureNode;

    double _fixedPressure;

    // for total area constraint
    GlobalConstraint _areaConstraint;
    double _area;
    double _prescribedArea;
    double _penaltyArea;
    MultiplierNode * _tensionNode;

    double _fixedTension;

    // for total curvature constraint
    GlobalConstraint _totalCurvatureConstraint;
    double _totalCurvature;
    double _prescribedTotalCurvature;
    double _penaltyTotalCurvature;
    MultiplierNode * _totalCurvatureNode;

    double _fixedTotalCurvatureForce;

    ConstraintContainer _constraints;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "C1AxiShellBody.icc"

#endif // __C1AxiShellBody_h__
