// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Capsid3DBody_h__)
#define __Capsid3DBody_h__

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
#include "NonlinearElastic.h"
#include "Contact.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  /*! 
    Capsid3DBody is a concrete class derived from Body, implementing
    the concept of a 3D body using tetrahedral finite elements.
  */
  template < class Quadrature_t, class Material_t, class Shape_t >
  class Capsid3DBody : public Body
  {
  public:
    
    // typedefs
    typedef 
    NonlinearElastic<Quadrature_t,Material_t,Shape_t> CapsidElement_t;
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
    std::vector< Constraint* > ContactContainer;
    typedef typename 
    ContactContainer::iterator ContactIterator;
    typedef typename 
    ContactContainer::const_iterator ConstContactIterator;
		
    typedef 
    std::vector<int> ElementConnectivity;
    typedef 
    std::vector<ElementConnectivity> ConnectivityContainer;
    typedef 
    ConnectivityContainer::iterator ConnectivityIterator;
    typedef 
    ConnectivityContainer::const_iterator ConstConnectivityIterator;


    //! Default Constructor
    Capsid3DBody() {;}

    //! Construct body from material, nodes, element connectivities, and quadrature order
    Capsid3DBody(
	Material_t material,
	ConnectivityContainer & connectivities,
	const NodeContainer & nodes,
	const unsigned quadOrder=1 ) 
    {
      initializeBody(material, connectivities, nodes, quadOrder);
    }

    //! Initialize body
    void initializeBody(
		Material_t material,
		ConnectivityContainer & connectivities,
		const NodeContainer & nodes,
		const unsigned quadOrder);
    
    //! Virtual destructor
    ~Capsid3DBody()
    {
      for(int i =0; i<_elements.size(); i++)
      {
	delete _elements[i];
      } 
    }
    
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

    //! loop over elements and add up the volume of each one
    double volume() const {
      double totalVolume = 0.0;
      for(ConstCapsidElementIterator e=_capsidElements.begin(); 
	  e!=_capsidElements.end(); e++) {
	totalVolume += (*e)->elVolume();
      }	
      return totalVolume;
    }
    
    //! loop over elements and add up the work of each one
    double work() const {
      double work = 0.0;
      for(ConstCapsidElementIterator e=_capsidElements.begin(); 
	  e!=_capsidElements.end(); e++) {
	work += -(*e)->work()+(*e)->constraintEnergy();
      }	
      return work;
    }

    //! Add a contact to body
    void pushBackContact( Constraint * c ) { _contacts.push_back( c ); }
    //! Print a Paraview file for a body made of linear tet elements
    void printParaviewLinearTet(const std::string name) const;
    //! Print a Paraview file for a body made of quadratic tet elements
    void printParaviewQuadTet(const std::string name) const;
    //! General printing of a Paraview file, calls either linear or quadratic print method based on number of nodes in an element
    void printParaview(const std::string name) const {
      int numNode = _capsidElements[0]->nodes().size();
      if(numNode == 4) printParaviewLinearTet(name);
      if(numNode == 10) printParaviewQuadTet(name);
      else std::cout << "Cannot print unless linear or quadratic." << std::endl;
    }
  
    //! Prints Paraview file with all stress/force information after post-processing
    void printParaviewPostProcess(const std::string name) const;
    //! For post-processing mechanics on body given a displacement field
    void PostProcess( const std::string name ) {
      compute( true, true, false );
      printParaviewPostProcess( name );      
    }

  private:

    //
    // data
    //
    
    //! Nodes
    CapsidNodeContainer _capsidNodes;

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

#include "Capsid3DBody.icc"

#endif // __Capsid3DBody_h__
