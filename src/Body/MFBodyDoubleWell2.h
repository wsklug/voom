// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Ankush Aggarwal & William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__MFBodyDoubleWell_h__)
#define __MFBodyDoubleWell_h__

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
#include "MFShape.h"
#include "Contact.h"
#include "voom.h"
#include "Node.h"
#include "NodeBin.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{
  
  /*! 
    MFBodyDoubleWell is a concrete class derived from Body, implementing
    the concept of a 3D body using mesh free RKP method and SNNI
  */
  template <class Material_t, class Shape_t >
  class MFBodyDoubleWell : public Body
  {
  public:
    
    // typedefs
    typedef DeformationNode<3> MFNode_t;
    typedef typename std::vector< MFNode_t* > MFNodeContainer;
    typedef typename MFNodeContainer::iterator MFNodeIterator;
    typedef typename MFNodeContainer::const_iterator ConstMFNodeIterator;
    
    typedef std::vector< Constraint* > ContactContainer;
    typedef typename ContactContainer::iterator ContactIterator;
    typedef typename ContactContainer::const_iterator ConstContactIterator;
    
    //! Default Constructor
    MFBodyDoubleWell() {;}
    
    //! Construct body from material, nodes and corresponding volumes
    MFBodyDoubleWell(
	   Material_t material,
	   const NodeContainer & nodes,
	   const std::vector<double> & node_volume,
	   const std::vector<double> & supp_size,
	   const std::vector<Tensor3D> & Fhatinv)
    {
      initializeBody(material, nodes, node_volume, supp_size, Fhatinv);
    }
    
    //! Construct body from material, nodes and corresponding volumes with modified kernel
    MFBodyDoubleWell(
	   Material_t material,
	   const NodeContainer & nodes,
	   const std::vector<double> & node_volume,
	   const std::vector<double> & supp_size,
	   const std::vector<double> & supp_size_hat, 
	   const std::vector<Tensor3D> & Fhatinv)
    {
      initializeBody(material, nodes, node_volume, supp_size,supp_size_hat,Fhatinv);
    }
    //! Initialize body with cubic-kernel
    void initializeBody(
			Material_t material,
			const NodeContainer & nodes,
			const std::vector<double> & node_volume,
			const std::vector<double> & supp_size,
			const std::vector<Tensor3D> & Fhatinv);
    
    //! Initialize body with modified kernel
    void initializeBody(
			Material_t material,
			const NodeContainer & nodes,
			const std::vector<double> & node_volume,
			const std::vector<double> & supp_size,
			const std::vector<double> & supp_size_hat,
			const std::vector<Tensor3D> & Fhatinv);
    
    //! Initialize body with singular kernel
    void initializeBody(
			Material_t material,
			const NodeContainer & nodes,
			const std::vector<double> & node_volume,
			const std::vector<double> & supp_size,
			const double sing_order,
			const std::vector<Tensor3D> & Fhatinv);
    //! Virtual destructor
    virtual ~MFBodyDoubleWell() {;}
    
    //! Do mechanics on Body
    void compute( bool f0, bool f1, bool f2 );
    
    //! Compute the invariants
    void cal_invariants(std::vector<double> & I1, std::vector<double> & I2, std::vector<double> & I3);
    
    //! Compute the stretch ratios along principal directions
    void cal_stretch_ratios( std::vector<double> & lam1, std::vector<double> & lam2, std::vector<double> & lam3 );
    //! Return the energy of the body
    double totalStrainEnergy() const { return _energy;}
    
    /*    //! loop over nodes and add up the strain energy of each one
	  double totalStrainEnergy() const {
	  double totalStrainEnergy = 0.0;
	  
	  for(QuadPointIterator p=_quadPoints.begin(); 
	  p!=_quadPoints.end(); p++) {
	  totalStrainEnergy += p->strainEnergy();
	  }	
	  return totalStrainEnergy;
	  }*/

    //! loop over nodes and add up the volume of each one
    double volume() const {
      double totalVolume = 0.0;
      for(ConstQuadPointIterator p=_quadPoints.begin(); 
	  p!=_quadPoints.end(); p++) {
	totalVolume += p->weight;
      }	
      return totalVolume;
    }
    
    //! loop over nodes and add up the work of each one
    double work() const {
      double work = 0.0;
      
      for(QuadPointIterator p=_quadPoints.begin(); 
	  p!=_quadPoints.end(); p++) {
	work += -p->work()+p->constraintEnergy();
      }	
      return work;
    }
    
    //! Add a contact to body
    void pushBackContact( Constraint * c ) { _contacts.push_back( c ); }
    
    //! General printing of a Paraview file
    void printParaview(const std::string name) const;// {printParaviewLinearTet(name);}
    
    //! Prints Paraview file with all stress/force information after post-processing
    void printParaviewPostProcess(const std::string name) const; //{};
    
    //! For post-processing mechanics on body given a displacement field
    void PostProcess( const std::string name ) {
      compute( true, true, false );
      printParaviewPostProcess( name );      
    }
    
    //quadrature points structurte
    struct QuadPointStruct 
    {
      double	 weight;
      Tensor3D Fhatinv;
      typename Shape_t::FunctionContainer shapeFunctions;
      typename Shape_t::FunctionContainer shapexDerivatives, shapeyDerivatives, shapezDerivatives;
      //neighbours list
      typename Shape_t::NodeNContainer neighbours;
      Material_t material;
      QuadPointStruct(double w, const Material_t & m, const Shape_t & s, Tensor3D F);
    };
    
    typedef typename std::vector<QuadPointStruct> QuadPointContainer;
    typedef typename QuadPointContainer::iterator QuadPointIterator;
    typedef typename QuadPointContainer::const_iterator ConstQuadPointIterator;
    
    //! Access the container of quadrature points
    const QuadPointContainer & quadraturePoints() const { return _quadPoints; }
    QuadPointContainer & quadraturePoints() { return _quadPoints; }
    
    
  private:
    
    //
    // data
    //
    
    //! Nodes
    MFNodeContainer _MFnodes;
    
    //! Contacts
    ContactContainer _contacts;
    
    //! Container of quadrature points
    QuadPointContainer _quadPoints;
    
#ifdef WITH_MPI
    int _processorRank;
    int _nProcessors;
#endif
    
  };  
} // namespace voom

#include "MFBodyDoubleWell.icc"

#endif // __MFBodyDoubleWell_h__
