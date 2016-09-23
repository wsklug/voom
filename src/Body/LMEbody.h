// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Ankush Aggarwal & William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__LMEbody_h__)
#define __LMEbody_h__

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
#include "LMEshape.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

// External LAPACK routine to invert matrices 
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                        double *b, int *ldb, int *info);

namespace voom
{
  
  /*! 
    LMEBody is a concrete class derived from Body, implementing
    the concept of a 3D body using mesh free LME method
  */
  template <class Material_t, class Shape_t >
  class LMEbody : public Body
  {
  public:
    
    // typedefs
    typedef DeformationNode<3> MFNode_t;
    typedef typename std::vector< MFNode_t* > NodeContainer;
    typedef typename NodeContainer::iterator NodeIterator;
    typedef typename NodeContainer::const_iterator ConstNodeIterator;
    
    //! Default Constructor
    LMEbody() {}
    
    //! Construct body from material, nodes, volumes, and LME parameters
    LMEbody(Material_t material,
	    const NodeContainer & nodes,
	    const std::vector<double> & node_volume,
	    const std::vector<double> & supp_size,
	    const double beta, const double searchR, const double tol, const int nItMax);

    
    //! Destructor
    ~LMEbody() 
    {
      for (QuadPointIterator p =_quadPoints.begin(); p!= _quadPoints.end(); p++)
      {
	free(p->A);
      }
    }
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! Compute the invariants
    void cal_invariants(std::vector<double> & I1, std::vector<double> & I2, std::vector<double> & I3);

    //! Return the energy of the body
    double totalStrainEnergy() const { return _energy;}

    //! Return the body volume 
    double volume() const 
    {
      double totalVolume = 0.0;
      for(ConstQuadPointIterator p = _quadPoints.begin(); 
	  p != _quadPoints.end(); p++) {
	totalVolume += p->weight;
      }	
      return totalVolume;
    }
    
    //! General printing of a Paraview file
    void printParaview(const std::string name) const;
    
   
    
    //quadrature points structurte
    struct QuadPointStruct 
    {
      double	 weight;
      typename Shape_t::FunctionContainer shapeFunctions;
      typename Shape_t::FunctionContainer shapexDerivatives, shapeyDerivatives, shapezDerivatives;
      // neighbours list
      typename Shape_t::NodeNContainer neighbours;
      Material_t material;
      double *A;
      QuadPointStruct(double w, const Material_t & m, const Shape_t & s, double *A);
    };
    
    typedef typename std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }
    QuadPointContainer & quadraturePoints() { return _quadPoints; }
    
    
  private:
    // data
    
    //! Nodes
    const NodeContainer & _nodes;
    
    //! Container of quadrature points
    QuadPointContainer _quadPoints;

    //! Body volume
    double _volume;

    //! Stored energy
    double _energy;

    //! Shape functions parameters
    const double _beta;
    const double _searchR;
    const double _tol;
    const int _nItMax; 
    
#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
}

#include "LMEbody.cc"

#endif // __LMEbody_h__
