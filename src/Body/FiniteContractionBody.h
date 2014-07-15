// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   Mainak Sarkar & William S. Klug
//                University of California Los Angeles
//                   (C) 2008 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__FiniteContractionBody_h__)
#define __FiniteContractionBody_h__

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
#include "FiniteContraction.h"
//#include "ads/halfedge.h"
#include "voom.h"
#include "Node.h"
//#include "HDSVertexNode.h"
//#include "HDSConnectivity.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    Capsid3DBody is a concrete class derived from Body, implementing
    the concept of a 3D capsid body using tetrahedral finite elements.
  */
template<class DefQuadrature_t, class Material_t, class DefShape_t,
	class VoltQuadrature_t, class EpMaterial_t, class VoltShape_t>
  class FiniteContractionBody : public Body
  {
  public:
    
    // typedefs
    typedef 
    FiniteContraction<DefQuadrature_t, Material_t, DefShape_t,
			VoltQuadrature_t, EpMaterial_t, VoltShape_t> ContractionElement_t;
    typedef typename 
    ContractionElement_t::Node_t ContractionNode_t;
    typedef typename 
    std::vector< ContractionNode_t* > ContractionNodeContainer;
    typedef typename 
    ContractionNodeContainer::iterator ContractionNodeIterator;
    typedef typename 
    ContractionNodeContainer::const_iterator ConstContractionNodeIterator;

    typedef typename 
    ContractionElement_t::Vnode_t ContractionVnode_t;
    typedef typename 
    std::vector< ContractionVnode_t* > ContractionVnodeContainer;
    typedef typename 
    ContractionVnodeContainer::iterator ContractionVnodeIterator;
    typedef typename 
    ContractionVnodeContainer::const_iterator ConstContractionVnodeIterator;

    typedef 
    std::vector< ContractionElement_t* > ContractionElementContainer;
    typedef typename 
    ContractionElementContainer::iterator ContractionElementIterator;
    typedef typename 
    ContractionElementContainer::const_iterator ConstContractionElementIterator;
		
    typedef 
    std::vector<int> ElementConnectivity;
    typedef 
    std::vector<ElementConnectivity> ConnectivityContainer;
    typedef 
    ConnectivityContainer::iterator ConnectivityIterator;
    typedef 
    ConnectivityContainer::const_iterator ConstConnectivityIterator;


    //! Default Constructor
    FiniteContractionBody() {;}

    //! Construct from stuff
    FiniteContractionBody(Material_t material, EpMaterial_t epMaterial,
			ConnectivityContainer & connectivities,
    		const ContractionNodeContainer & nodes,
    		const ContractionVnodeContainer & vNodes,
    		const unsigned int quadOrder)
    {
      initializeBody(material, epMaterial,
      		connectivities,nodes,vNodes, quadOrder);
    }

    //! initialize
    void initializeBody(Material_t material, EpMaterial_t epMaterial,
			ConnectivityContainer & connectivities,
			const ContractionNodeContainer & nodes,
			const ContractionVnodeContainer & vNodes,
			const unsigned int quadOrder);
    
    //! virtual destructor
    virtual ~FiniteContractionBody() {;}
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstContractionElementIterator e=_contractionElements.begin(); 
	  e!=_contractionElements.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
      }	
      return totalStrainEnergy;
    }
    
    //! loop over elements and add up the work of each one
    double work() const {
      double work = 0.0;
      for(ConstContractionElementIterator e=_contractionElements.begin(); 
	  e!=_contractionElements.end(); e++) {
	work += -(*e)->work()+(*e)->constraintEnergy();
      }	
      return work;
    }


    void printParaview(const std::string name) const{};
//    void printParaviewQuadTet(const std::string name) const;

  private:

    //! Nodes
    ContractionNodeContainer     _contractionNodes;
    // voltage nodes
    ContractionVnodeContainer     _contractionVnodes;
    //! Elements
    ContractionElementContainer  _contractionElements;		
    

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "FiniteContractionBody.icc"

#endif // __FiniteContractionBody_h___
