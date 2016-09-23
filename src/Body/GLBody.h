// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug, Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2007 All Rights Reserved
//
//----------------------------------------------------------------------



#if !defined(__GLBody_h__)
#define __GLBody_h__

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
#include "GLElement.h"
#include "Contact.h"
//#include "ads/halfedge.h"
#include "voom.h"
#include "Node.h"
#include "HDSVertexNode.h"
#include "HDSConnectivity.h"
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
  class GLBody : public Body
  {
  public:
    
    

    // typedefs
    typedef GLElement<Material_t> FeElement_t;
    typedef typename FeElement_t::Node_t FeNode_t;
    typedef typename std::vector< FeNode_t* > FeNodeContainer;
    typedef typename FeNodeContainer::iterator FeNodeIterator;
    typedef typename FeNodeContainer::const_iterator ConstFeNodeIterator;

    typedef std::vector< FeElement_t* > FeElementContainer;
    typedef typename FeElementContainer::iterator FeElementIterator;
    typedef typename FeElementContainer::const_iterator ConstFeElementIterator;
		
    typedef std::vector< Contact* > ContactContainer;
    typedef typename ContactContainer::iterator ContactIterator;
    typedef typename ContactContainer::const_iterator ConstContactIterator;

    typedef tvmet::Vector<int,3> ElementConnectivity;
    typedef std::vector<ElementConnectivity> ConnectivityContainer;

    //! Default Constructor
    GLBody() {;}

    //! Construct from stuff
    GLBody(Material_t material,
	   ConnectivityContainer & connectivities,
	   const NodeContainer & nodes,
	   const unsigned quadOrder = 1,
	   const int numberOfBoundaries = 0,
	   const double pressure = 0.0,
	   const double tension = 0.0,
	   const double penaltyVolume = 1.0e4,
	   const double penaltyArea= 1.0e6,
	   GlobalConstraint volumeConstraint = noConstraint,
	   GlobalConstraint areaConstraint = noConstraint ) {

      initializeBody(material, connectivities, nodes, quadOrder, 
		     numberOfBoundaries, pressure, tension,
		     penaltyVolume, penaltyArea,
		     volumeConstraint, areaConstraint);
    }

    //! initialize
    void initializeBody(Material_t material,
			ConnectivityContainer & connectivities,
			const NodeContainer &  nodes,
			const unsigned quadOrder,
			const int numberOfBoundaries,
			const double pressure,
			const double tension,
			const double penaltyVolume,
			const double penaltyArea,
			GlobalConstraint volumeConstraint,
			GlobalConstraint areaConstraint);

    
    //! virtual destructor
    virtual ~GLBody() {
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
// 	(*e)->setConstraintVolume( volumeFactor*(*e)->constraintVolume() );
// 	_prescribedVolume += (*e)->constraintVolume();
//       }

    }

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
	
      }	
      return totalStrainEnergy;
    }

    double bendingEnergy() const {
      double bendingEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	bendingEnergy += (*e)->bendingEnergy();
	
      }	
      return bendingEnergy;
    }

    double inplaneEnergy() const {
      double inplaneEnergy = 0.0;
      for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
	inplaneEnergy += (*e)->inplaneEnergy();
	
      }	
      return inplaneEnergy;
    }    
    //! loop over elements and add up the work of each one
    //double work() const {
    //double work = 0.0;
    //for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
    //work += -(*e)->work()+(*e)->constraintEnergy();
    //}	
    //return work;
    //}

//     void setCytoSpring(double new_mu, double new_kS, double new_kSpring){
//       for(ConstFeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
// 	(*e)->setCytoSpring(new_mu, new_kS, new_kSpring);
//       }	      
//     }


    //! print body connectivities by using Half-edge Data Structure info.
    void printByHDS();
        
    //! output to openDX format
    void createOpenDXData(const std::string& fileName, 
			  const int nWhichForce = 2);

    void printParaview(const std::string fileName) const ;

    void printObj(const std::string name) const;
    
    virtual void incrementLoads(const double load){}

    void resetReference() {
      for(FeNodeIterator n = _shellNodes.begin(); n!= _shellNodes.end(); n++) {
	(*n)->resetPosition();
      }

//       for(FeElementIterator e=_shells.begin(); e!=_shells.end(); e++) {
// 	(*e)->updateRefConfiguration();
//       }	
//      cytoskeleton energy not updated

      compute(true,true,true);//add this because we have nodal force add up twice (updateRefConfiguration() has compute)
                              //zero out force first. 
      return;
    }

    void pushBackContact( Contact * c ) { _contacts.push_back( c ); }

    //
    // Compute force along the generatix in the vescile (Paul's model)
    // for comparison
//     void ComputeVescileForce();
    
  private:

    //! Store corner valences in a vector for building irregular elements
    typedef tvmet::Vector<unsigned int, 3> CornerValences;
        
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

    // for total area constraint
    GlobalConstraint _areaConstraint;
    double _area;
    double _prescribedArea;
    double _penaltyArea;
    MultiplierNode * _tensionNode;

    ContactContainer _contacts;

    double _fixedPressure, _fixedTension;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "GLBody.cc"

#endif // __GLBody_h__
