// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__Body3D_h__)
#define __Body3D_h__

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
#include "Element3D.h"
#include "Quadrature.h"
#include "Shape.h"
#include "Contact.h"
#include "voom.h"
#include "Node.h"
#include "HomogMP.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{

  /*! 
    Body3D is a concrete class derived from Body, implementing
    the concept of a 3D body using tetrahedral finite elements 
    and computing energy invariants
  */
  
  class Body3D : public Body
  {
  public:

    typedef vector<DeformationNode<3> * > DefNodeContainer;
    typedef DefNodeContainer::iterator DefNodeContainerIterator;
    typedef DefNodeContainer::const_iterator ConstDefNodeContainerIterator;
		
    typedef vector<int> ElementConnectivity;
    typedef vector<ElementConnectivity> ConnectivityContainer;
    typedef ConnectivityContainer::iterator ConnectivityIterator;
    typedef ConnectivityContainer::const_iterator ConstConnectivityIterator;

    //! Default Constructor
    Body3D() {;}

    //! Construct body from material, element connectivities, nodes, quadrature and shape functions
    //! Initialize body
    Body3D(vector<Material *> Mat,
	   const ConnectivityContainer & Connectivities,
	   const DefNodeContainer & DefNodes,
	   Quadrature<3> & Quad,
	   Shape<3> & Sh,
	   double k = -1.0);
    
    //! Virtual destructor
    ~Body3D()
    {
      for(int i =0; i<_elements.size(); i++)
      {
	delete _elements[i];
      } 
    }
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! loop over elements and add up the strain energy of each one
    double totalStrainEnergy() const 
    {
      double totalStrainEnergy = 0.0;
      for(ConstElementIterator e = _elements.begin(); e!=_elements.end(); e++)
      {
	// assert(dynamic_cast<Element3D* >(*e) != NULL);
	totalStrainEnergy += dynamic_cast<Element3D* >(*e)->strainEnergy();
      }	

      return totalStrainEnergy;
    }

    //! loop over elements and add up the volume of each one
    double volume() const
    {
      double totalVolume = 0.0;
      for(ConstElementIterator e = _elements.begin(); e != _elements.end(); e++)
      {
	totalVolume += dynamic_cast<Element3D* >(*e)->volume();
      }	

      return totalVolume;
    }

    vector<pair<Vector3D, vector<double > > > invariants(int &f)
    {
      vector<pair<Vector3D, vector<double > > > Invariants(f);
      for(ConstElementIterator e = _elements.begin(); e != _elements.end(); e++)
      {
	vector<pair<Vector3D, vector<double > > > invariantsLoc = dynamic_cast<Element3D* >(*e)->invariants(f);
	Invariants.insert(Invariants.end(), invariantsLoc.begin(), invariantsLoc.end()); 
      }	

      return Invariants;
    };

    int invariants(vector<double > &I1, vector<double > &I2, vector<double > &I3)
    {
      int f = 0;
      // Works only with 1 QP per element !!
      for(int el = 0; el < _elements.size(); el++)
      {
	vector<pair<Vector3D, vector<double > > > invariantsLoc = dynamic_cast<Element3D* >(_elements[el])->invariants(f);
	vector<double > Inv = invariantsLoc[0].second;
	I1[el] = Inv[0];
	I2[el] = Inv[1];
	I3[el] = Inv[2];
      }	
      return f;

    };
    
    //! set I1mostProb (most probable invariants) etc. for the experimental material which helps to find the "true" reference state
    void setMPinv(vector<double> & I1mp, vector<double> & I2mp, vector<double> & Jmp);

    //! set I1mostProb (most probable invariants) etc. for the experimental material which helps to find the "true" reference state
    void reset();



    //! Print a Paraview file for a body made of linear tet elements
    void printParaviewLinearTet(const string name) const;
    //! Print a Paraview file for a body made of quadratic tet elements
    void printParaviewQuadTet(const string name) const;
    //! General printing of a Paraview file, calls either linear or quadratic print method based on number of nodes in an element
    void printParaview(const string name) const
    {/*
      int numNode = dynamic_cast<Element3D* >(_elements[0])->nodes().size();
      if(numNode == 4) printParaviewLinearTet(name);
      if(numNode == 10) printParaviewQuadTet(name);
      else cout << "Cannot print unless linear or quadratic." << endl;*/
    }
  
    //! Prints Paraview file with all stress/force information after post-processing
    void printParaviewPostProcess(const string name) const;
    //! For post-processing mechanics on body given a displacement field
    void PostProcess( const string name )
    {
      compute( true, true, false );
      printParaviewPostProcess( name );      
    }

  private:
    //
    // data
    //
   

#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "Body3D.cc"

#endif // __Body3D_h__
