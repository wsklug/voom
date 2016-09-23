// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug
//                University of California Los Angeles
//                   (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__CapsidAtomisticBody_h__)
#define __CapsidAtomisticBody_h__

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
#include "AtomisticCoupled.h"
#include "Contact.h"
#include "ads/halfedge.h"
#include "voom.h"
#include "Node.h"
#include "HDSVertexNode.h"
#include "HDSConnectivity.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    CapsidAtomisticBody is a concrete class derived from Body, implementing
    the concept of a capsid body using hexahedral finite elements, and 
    calculating forces at the atomic level.
  */
  template < class Shape_t >
  class CapsidAtomisticBody : public Body
  {
  public:
    
    // typedefs
    typedef 
      AtomisticCoupled<Shape_t> CapsidElement_t;
    typedef typename 
      CapsidElement_t::Node_t CapsidNode_t;
    typedef typename 
      std::vector< CapsidNode_t* > CapsidNodeContainer;
    typedef typename 
      CapsidNodeContainer::iterator CapsidNodeIterator;
    typedef typename 
      CapsidNodeContainer::const_iterator ConstCapsidNodeIterator;
    
    typedef 
      std::vector< CapsidElement_t* > CapsidElementContainer;
    typedef typename 
      CapsidElementContainer::iterator CapsidElementIterator;
    typedef typename 
      CapsidElementContainer::const_iterator ConstCapsidElementIterator;
    
    typedef 
      std::vector< Contact* > ContactContainer;
    typedef typename 
      ContactContainer::iterator ContactIterator;
    typedef typename 
      ContactContainer::const_iterator ConstContactIterator;
    
    typedef 
    std::vector< DeformationNode<3>* > CapsidAtomContainer;
    typedef
      CapsidAtomContainer::iterator CapsidAtomIterator;
    typedef
      CapsidAtomContainer::const_iterator ConstCapsidAtomIterator;
           
    typedef 
      std::vector<int> ElementConnectivity;
    typedef 
      std::vector<ElementConnectivity> ConnectivityContainer;
    typedef 
      ConnectivityContainer::iterator ConnectivityIterator;
    typedef 
      ConnectivityContainer::const_iterator ConstConnectivityIterator;


    //! Default Constructor
    CapsidAtomisticBody() {;}

    //! Construct from stuff
    CapsidAtomisticBody(
	ConnectivityContainer & connectivities,
	const NodeContainer & nodes,
	std::vector<int> & atomElementConnectivity,
	CapsidAtomContainer & atoms) 
    {
      initializeBody(connectivities, 
		     nodes, 
		     atomElementConnectivity, 
		     atoms);
    }

    //! initialize
    void initializeBody(
		ConnectivityContainer & connectivities,
		const NodeContainer & nodes,
		std::vector<int> & atomElementConnectivity,
		CapsidAtomContainer & atoms);
    
    //! virtual destructor
    virtual ~CapsidAtomisticBody() {;}
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const {
      double totalStrainEnergy = 0.0;
      for(ConstCapsidElementIterator e=_capsidElements.begin(); 
	  e!=_capsidElements.end(); e++) {
	totalStrainEnergy += (*e)->strainEnergy();
      }	
      return totalStrainEnergy;
    }

    void pushBackContact( Contact * c ) { _contacts.push_back( c ); }

    void printParaview(const std::string name) const;
    // need to make this an empty function since it's virtual in Body.h
    void printParaviewQuadTet(const std::string name) const {;}
    
/*     // For post-processing mechanics on body given a displacement field */
/*     void printParaviewPostProcess(const std::string name) const; */
/*     void PostProcess( const std::string name ) { */
/*       compute( true, true, false ); */
/*       printParaviewPostProcess( name );       */
/*     } */

    CapsidAtomContainer & capsidAtoms() {return _capsidAtoms;}
    CapsidElementContainer & capsidElements() {return _capsidElements;}

  private:

    //
    // data
    //
    
    //! Nodes
    CapsidNodeContainer _capsidNodes;

    //! Atoms
    CapsidAtomContainer _capsidAtoms;

    //! Elements
    CapsidElementContainer 	_capsidElements;		

    //! Contacts
    ContactContainer _contacts;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "CapsidAtomisticBody.icc"

#endif // __CapsidAtomisticBody_h__
