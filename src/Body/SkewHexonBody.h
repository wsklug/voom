// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Ankush Aggarwal
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__SkewHexonBody_h__)
#define __SkewHexonBody_h__

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
#include "ModEvansElastic.h"
#include "EvansElastic_Skewed.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    SkewHexonBody is a concrete class derived from Body similar to C0MembraneBody,
    implementing the concept of a membrane body using C0-conforming finite elements
    with EvansElastic_Skewed having skewed equilibrium state.
  */
  template < class Quadrature_t, class Shape_t >
  class SkewHexonBody : public Body
  {
  public:
    
    // typedefs
    typedef 
    C0Membrane<Quadrature_t,EvansElastic_Skewed,Shape_t> MembraneElement_t;
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
    SkewHexonBody() {;}

    //! Construct from stuff
    SkewHexonBody(
	//Material_t material instead pass material properties
	ConnectivityContainer & connectivities,
	const NodeContainer & nodes,
        const std::vector<double> & shear_angle,
        const double eta,
        const double mu,
        const double kS,
        const double c1,
        const double c2,
	const unsigned quadOrder=1,
	const double pressure=0.0, 
	const double tension=0.0,
	const double penaltyVolume=1.0e4,
	const double penaltyArea=1.0e4,
	GlobalConstraint volumeConstraint = noConstraint,
	GlobalConstraint areaConstraint = noConstraint) ;
    
    //! virtual destructor
    virtual ~SkewHexonBody() {;}
    
//     void setMaterialY(double Y) {
//       for(int si=0; si<_membranes.size(); si++) {
// 	MembraneElement_t* s=_membranes[si];
// 	s->setMaterialY(Y);     	
//       }
//     }

    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! Set the reference configuration to an equilateral triangle in XY plane    
    void SetRefConfiguration(double edge);

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

    void printParaview(const std::string fileName) const ;

    void printParaview2(const std::string fileName) const ;

    void printParaviewEigVec(const std::string name) const;

    void printObj(const std::string fileName) const ;
    
    //! recreate input file
    void createInputFile(const std::string& filename);
    
    void resetReference() {
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

    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    //! Query an element's activity status
    bool active(int e) { return _active[e]; }

    //! Mark an element as active so it will be computed
    void activate(int e) { _active[e] = true; }

    //! Mark an element as inactive so it will not be computed
    void deactivate(int e) { _active[e] = false; }

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

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "SkewHexonBody.icc"

#endif // __SkewHexonBody_h__
