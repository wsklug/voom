// -*- C++ -*-
//----------------------------------------------------------------------
//
//           Ankush Aggarwal, William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__StretchHexonBody_h__)
#define __StretchHexonBody_h__

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
#include "C0MembraneStretch.h"
#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include "EvansElastic_Stretch.h"
#include "ShapeTri3.h"
#include "Quadrature.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{

  /*! 
    StretchHexonBody is a concrete class derived from Body similar to C0MembraneBody,
    implementing the concept of a membrane body using C0-conforming finite elements
    with EvansElastic_Stretch having stretch equilibrium state.
  */
  class StretchHexonBody : public Body
  {
  public:
    
    // typedefs
    typedef DeformationNode<3> MembraneNode;
    typedef vector< MembraneNode* > MembraneNodeContainer;
    typedef MembraneNodeContainer::iterator MembraneNodeIterator;
    typedef MembraneNodeContainer::const_iterator ConstMembraneNodeIterator;

    typedef vector< ScalarFieldNode<3>* > StretchNodeContainer;
    typedef StretchNodeContainer::iterator StretchNodeIterator;
    typedef StretchNodeContainer::const_iterator ConstStretchNodeIterator;

    typedef C0MembraneStretch MembraneElement;
    typedef vector<MembraneElement* > MembraneElementContainer;
    typedef MembraneElementContainer::iterator MembraneElementIterator;
    typedef MembraneElementContainer::const_iterator ConstMembraneElementIterator;

    typedef vector<int> ElementConnectivity;
    typedef vector<ElementConnectivity> ConnectivityContainer;
    typedef ConnectivityContainer::iterator ConnectivityIterator;
    typedef ConnectivityContainer::const_iterator ConstConnectivityIterator;
		
    typedef vector<Constraint* > ConstraintContainer;
    typedef ConstraintContainer::iterator ConstraintIterator;
    typedef ConstraintContainer::const_iterator ConstConstraintIterator;
		
   


    //! Default Constructor
    StretchHexonBody() {;}

    //! Constructor
    StretchHexonBody(const ConnectivityContainer & connectivities,
		     const MembraneNodeContainer & DefNodes,
		     const StretchNodeContainer & StretchNodes,
		     const StretchNodeContainer & DirectionNodes,
		     const vector<double> & stretch_angle,
		     const double mu,
		     const double kS,
		     const int WcType,
		     double* WcConst,
		     const Quadrature<2> & Quad,
		     Shape<2> * shape,
		     const vector<unsigned int> & CapsomersNum,
		     ScalarFieldNode<3>* phiNode = NULL);
    
    //! virtual destructor
    ~StretchHexonBody()
    {
      for (MembraneElementIterator mel = _membraneElements.begin();  mel != _membraneElements.end(); mel++)
      {
	MembraneElement::QuadPointContainer quad = (*mel)->quadPoints();
	for (MembraneElement::QuadPointIterator q = quad.begin();  q != quad.end(); q++)
	{
	  delete (q->material);	  
	}
	delete (*mel);
      } 
    }
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! Set the reference configuration to an equilateral triangle in XY plane    
    void SetRefConfiguration(double edge)
    {
      for(int i = 0 ; i < _membraneElements.size(); i++) {
	if( !_active[i] ) continue;
	_membraneElements[i]->SetRefConfiguration(edge);
      }
    }

    void resetReference() {
      for(MembraneNodeIterator n = _membraneNodes.begin(); n != _membraneNodes.end(); n++)
	(*n)->resetPosition();

      for(MembraneElementIterator e = _membraneElements.begin();  e != _membraneElements.end(); e++)
	(*e)->updateRefConfiguration();
      
      return;
    }



    //! reset angle offset in all elements
    void setAngleOffset(vector<double > AngleOffset)
    {
      for (unsigned int i=0; i<_membraneElements.size(); i++)
	_membraneElements[i]->setAngleOffset(AngleOffset[i]);
    }
    


    //! get total area and volume
    double area() const { return _area; }
    double volume() const{ return _volume; }

    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const
    {
      double totalStrainEnergy = 0.0;
      for(ConstMembraneElementIterator e = _membraneElements.begin(); e != _membraneElements.end(); e++)
      {
	totalStrainEnergy += (*e)->strainEnergy();
	// cout << (*e)->strainEnergy() << endl;
      }	
      return totalStrainEnergy;
    }

    double totalConformationalEnergy() const
    {
      double totalConformationalEnergy = 0.0;
      for(ConstMembraneElementIterator e = _membraneElements.begin(); e != _membraneElements.end(); e++)
      {
	totalConformationalEnergy += (*e)->conformationalEnergy();
	// cout << (*e)->strainEnergy() << endl;
      }	
      return totalConformationalEnergy;
    }
    
    virtual void pushBackMemebranbeElement( MembraneElement* e ) { 
      _elements.push_back(e); 
      _membraneElements.push_back(e);
      _active.push_back(true);
    }



    void printParaview(const std::string fileName) const ;

    void printParaview2(const std::string fileName) const ;

    void printParaviewEigVec(const std::string name) const;
    
    //! recreate input file
    void createInputFile(const std::string& filename);



    void pushBackConstraint( Constraint * c ) { _constraints.push_back( c ); }

    //! Query an element's activity status
    bool active(int e) { return _active[e]; }

    //! Mark an element as active so it will be computed
    void activate(int e) { _active[e] = true; }

    //! Mark an element as inactive so it will not be computed
    void deactivate(int e) { _active[e] = false; }

  private:
    // data

    //! Nodes
    vector<DeformationNode<3>* >  _membraneNodes;
    StretchNodeContainer _stretchNodes;
    StretchNodeContainer _directionNodes;
      // ScalarFieldNode<3>* _phiNode; 

    //! Elements
    MembraneElementContainer _membraneElements;
    //! vector of boolean flags to activate/deactivate elements
    std::vector<bool> _active;

    double _area;
    double _volume;

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "StretchHexonBody.cc"

#endif // __StretchHexonBody_h__
