// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__C0MembraneBody_h__)
#define __C0MembraneBody_h__

#include <blitz/array.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Body.h"
#include "C0Membrane.h"
#include "Constraint.h"
#include "HalfEdgeMesh.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    C0MembraneBody is a concrete class derived from Body, implementing
    the concept of a membrane body using C0-conforming finite elements.
  */
  template < class Quadrature_t, class Material_t, class Shape_t >
  class C0MembraneBody : public Body
  {
  public:
    
    // typedefs
    typedef 
    C0Membrane<Quadrature_t,Material_t,Shape_t> MembraneElement_t;
    typedef typename 
    MembraneElement_t::Node_t MembraneNode_t;
    typedef typename 
    std::vector< MembraneNode_t* > MembraneNodeContainer;
    typedef typename 
    MembraneNodeContainer::iterator MembraneNodeIterator;
    typedef typename 
    MembraneNodeContainer::const_iterator ConstMembraneNodeIterator;

    typedef 
    std::vector< MembraneElement_t* > MembraneElementContainer;
    typedef typename 
    MembraneElementContainer::iterator MembraneElementIterator;
    typedef typename 
    MembraneElementContainer::const_iterator ConstMembraneElementIterator;
		
    typedef 
    std::vector< Constraint* > ConstraintContainer;
    typedef typename 
    ConstraintContainer::iterator ConstraintIterator;
    typedef typename 
    ConstraintContainer::const_iterator ConstConstraintIterator;
		
    typedef 
    std::vector<int> ElementConnectivity;
    typedef 
    std::vector<ElementConnectivity> ConnectivityContainer;
    typedef 
    ConnectivityContainer::iterator ConnectivityIterator;
    typedef 
    ConnectivityContainer::const_iterator ConstConnectivityIterator;


    //! Default Constructor
    C0MembraneBody() {;}

    //! Construct from stuff
    C0MembraneBody(
	Material_t material,
	ConnectivityContainer & connectivities,
	const NodeContainer & nodes,
	const unsigned quadOrder=1,
	const double pressure=0.0, 
	const double tension=0.0,
	const double penaltyVolume=1.0e4,
	const double penaltyArea=1.0e4,
	GlobalConstraint volumeConstraint = noConstraint,
	GlobalConstraint areaConstraint = noConstraint) ;
    
    //! virtual destructor
    virtual ~C0MembraneBody() {;}
    
//     void setMaterialY(double Y) {
//       for(int si=0; si<_membranes.size(); si++) {
// 	MembraneElement_t* s=_membranes[si];
// 	s->setMaterialY(Y);     	
//       }
//     }

    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    double volume() const{ return _volume; }
    
    //! get total area
    double area() const { return _area; }

    //! get constraint area
    double prescribedArea() const { return _prescribedArea; }

    //! get constraint volume
    double prescribedVolume() const { return _prescribedVolume; }
    
    //! get pressure
    double pressure() const {return _pressureNode->point();}

    //! set pressure
    void setPressure(double p) { _pressureNode->setPoint(p);}

    //! set constraint area
    void setConstraintArea(double A) { _prescribedArea = A;}

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstMembraneElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
      }	
      return totalStrainEnergy;
    }
    
    //! loop over elements and add up the work of each one
    double work() const {
      double work = 0.0;
      for(ConstMembraneElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	work += -(*e)->work()+(*e)->constraintEnergy();
      }	
      return work;
    }

    virtual void pushBack( Element* e ) { 
      _elements.push_back(e); 
      _active.push_back(true);}

    void updateFixedTension( double T ){
      _fixedTension  = T;
    }

    void updatePenaltyArea( double pA ){
      _penaltyArea = pA;
    }

    void printParaview(std::string fileName) const ;

    void printObj(const std::string fileName) const ;
    
    //! recreate input file
    void createInputFile(const std::string& filename);
    
    void resetReference() {
      // reset node ref config to current config and update elements
      // from there
      for(MembraneNodeIterator n = _membraneNodes.begin(); 
	  n != _membraneNodes.end(); n++) {
	(*n)->resetPosition();
      }
      
      for(MembraneElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	(*e)->updateRefConfiguration();
      }	
      return;
    }

    //! Reset element ref configs to equilateral with average edge length.
    void resetEquilateral() {
      // compute average element edge length
      double hAve = 0.0;
      for(int e = 0; e<_membranes.size(); e++) {
	double h=0.0;
	const MembraneNodeContainer & nds = _membranes[e]->nodes();
	for(int a=0; a<nds.size(); a++) {
	  h += norm2( nds[a]->position() - nds[(a+1) % nds.size()]->position() );
	}
	h /= nds.size();
	hAve += h;
      }
      hAve /= _membranes.size();

      std::cout << "Reseting all elements to be equilateral with element edge length " << hAve << std::endl; 
      
      // reset all elements to be equilalateral with average edge
      // length in ref config.
      for(MembraneElementIterator e=_membranes.begin(); 
	  e!=_membranes.end(); e++) {
	(*e)->updateRefConfiguration( hAve );
      }	

      return;
    }

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    //! Query an element's activity status
    bool active(int e) { return _active[e]; }

    //! Mark an element as active so it will be computed
    void activate(int e) { _active[e] = true; }

    //! Mark an element as inactive so it will not be computed
    void deactivate(int e) { _active[e] = false; }
    
    void checkElementConsistency() {
      for(int e=0; e<_membranes.size(); e++) {
	_membranes[e]->checkConsistency();
      }
    }

    void addStressDirection(const Vector3D & N) {
      _stressDirections.push_back(N);
    }

  private:

    //
    // data
    //
    
    //! Nodes
    MembraneNodeContainer _membraneNodes;

    //! Elements
    MembraneElementContainer 	_membranes;		

    //! vector of boolean flags to activate/deactivate elements
    std::vector<bool> _active;

    GlobalConstraint _volumeConstraint;
    GlobalConstraint _areaConstraint;

    double _penaltyArea;
    double _penaltyVolume;

    double _fixedPressure;
    double _fixedTension;

    double _area;
    double _prescribedArea;
    MultiplierNode * _tensionNode;

    double _volume;
    double _prescribedVolume;
    MultiplierNode * _pressureNode;

    ConstraintContainer _constraints;

    //! Reference Normal vectors to compute stresses along
    std::vector<Vector3D> _stressDirections;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "C0MembraneBody.icc"

#endif // __C0MembraneBody_h__
