// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2008 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file SemiflexibleGel.h

  \brief SemiflexibleGel is a concrete class derived from Body, implementing
  the concept of a collection of cross-linked semiflexible polymers (i.e., beams)

*/

#if !defined(__nlSemiflexibleGel_h__)
#define __nlSemiflexibleGel_h__

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

#include "EntropicSpring.h"
#include "AngleSpring.h"
#include "BrownianRod.h"

#include "Constraint.h"
#include "voom.h"
#include "Node.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Crosslink.h"
#include "Motor.h"
#include "IntersectionFinder.h"
#include "PeriodicTie.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{


  /*!  
    Concrete class for a semiflexible gel
  */
  template<int N>
  class SemiflexibleGel : public Body
  {
  public:
    
    // typedefs
    typedef BrownianNode<N> DefNode;
    typedef typename std::vector< DefNode* > DefNodeContainer;
    typedef typename DefNodeContainer::iterator DefNodeIterator;
    typedef typename DefNodeContainer::const_iterator ConstDefNodeIterator;

    typedef EntropicSpring<N> Bond;
    typedef std::vector< Bond* > BondContainer;
    typedef typename BondContainer::iterator BondIterator;
    typedef typename BondContainer::const_iterator ConstBondIterator;
		
    typedef AngleSpring<N> Angle;
    typedef std::vector< Angle* > AngleContainer;
    typedef typename AngleContainer::iterator AngleIterator;
    typedef typename AngleContainer::const_iterator ConstAngleIterator;

    typedef BrownianRod<N> Rod;
    typedef std::vector< Rod* > RodContainer;
    typedef typename RodContainer::iterator RodIterator;
    typedef typename RodContainer::const_iterator ConstRodIterator;
		
    typedef std::vector< Constraint* > ConstraintContainer;
    typedef typename ConstraintContainer::iterator ConstraintIterator;
    typedef typename ConstraintContainer::const_iterator ConstConstraintIterator;
		
    struct Filament {
      
      Filament( const DefNodeContainer & n, double kAngle, double viscosity, 
		double kT, double dt, double kC, int fitOrder);

      ~Filament();

      DefNodeContainer nodes;
      BondContainer    bonds;
      AngleContainer   angles;
      RodContainer     rods;

    };

    typedef std::vector< Filament* > FilamentContainer;
    typedef typename FilamentContainer::iterator FilamentIterator;
    typedef typename FilamentContainer::const_iterator ConstFilamentIterator;
    
    //Mo	
    typedef Crosslink<N> Clink;
    typedef std::vector< Clink* > CrosslinkContainer;
    typedef typename CrosslinkContainer::iterator CrosslinkIterator;
    typedef typename CrosslinkContainer::const_iterator ConstCrosslinkIterator;
  
    typedef Motor<N> MolMot;
    typedef std::vector< MolMot* > MotorContainer;
    typedef typename MotorContainer::iterator MotorIterator;
    typedef typename MotorContainer::const_iterator ConstMotorIterator;
		
    typedef tvmet::Vector<double,N> VectorND;

    typedef typename std::map< DefNode*, DefNode* > CrosslinkNodeMap;

    //! Default Constructor
    SemiflexibleGel() {
      _output=paraview;
    }

    //! Another Constructor
    SemiflexibleGel(DefNodeContainer & dNodes, VectorND & size, double filDens, int nodesPerFil, double nodeLen, double kAngle, double visc, double kT, double dt, double kcl, double shear, double kC, int fitOrder) {
      setRNGNorm(0.0,sqrt(kT/kAngle));
      setRNGUni();
      double vol=1.0;
      for(int i=0; i<N; i++) {
	vol *= size[i];
      }
      if(N == 2) {
	if(shear == 0.0) {
	  _box = new PeriodicBox(size[0],size[1]); 
	  std::cout << "Using simple periodic boundary conditions." << std::endl;
	}	  
	else {
	  _box = new LeesEdwards(size[0],size[1],shear);
	  std::cout << "Using Lees-Edwards boundary conditions." << std::endl;
	}
      }
      int nFils = (int)(filDens*vol);
      int id;
      NodeBase::DofIndexMap idx(N);
      DefNode * newDN;
      VectorND startPos;
      VectorND oldVec;
      VectorND newVec;
      double newAng;
      for(int i=0; i<nFils; i++) {
	_rnguni->seed((unsigned int)time(0)+i);
	_rngnorm->seed((unsigned int)time(0)+i);
	DefNodeContainer tmpDNC(nodesPerFil);
	// place the first node in the filament //
	id = i*nodesPerFil;
	for(int k=0; k<N; k++) {
	  idx[k] = N*id + k;
	  startPos[k] = size[k]*((_rnguni)->random());
	}
	newDN = new BrownianNode<N>(id,idx,startPos,startPos);
	newDN->setId(id);
	tmpDNC[0] = newDN;
	dNodes.push_back(newDN);
	// place the second node at some random angle //
	id++;
	for(int k=0; k<N; k++) {
	  idx[k] = N*id + k;
	}
	newAng = 2.0*3.14159*((_rnguni)->random());
	startPos[0] = nodeLen*cos(newAng) + (tmpDNC[0]->point())[0];
	startPos[1] = nodeLen*sin(newAng) + (tmpDNC[0]->point())[1];
	newDN = new BrownianNode<N>(id,idx,startPos,startPos);
	newDN->setId(id);
	tmpDNC[1] = newDN;
	dNodes.push_back(newDN);
	// now place the rest of the nodes, choosing the angles from a Boltzmann distribution //
	for(int j=2; j<nodesPerFil; j++) {
	  _rnguni->seed((unsigned int)time(0)+i*nodesPerFil+j);
	  _rngnorm->seed((unsigned int)time(0)+i*nodesPerFil+j);
	  id++;
	  for(int k=0; k<N; k++) {
	    idx[k] = N*id + k;
	  }
	  newAng = ((_rngnorm)->random());
	  oldVec = tmpDNC[j-1]->point() - tmpDNC[j-2]->point();
	  newVec[0] = oldVec[0]*cos(newAng) + oldVec[1]*sin(newAng);
	  newVec[1] = oldVec[1]*cos(newAng) - oldVec[0]*sin(newAng);
	  startPos = tmpDNC[j-1]->point() + newVec;
	  newDN = new BrownianNode<N>(id,idx,startPos,startPos);
	  newDN->setId(id);
	  tmpDNC[j] = newDN;
	  dNodes.push_back(newDN);
	}
	// add the filament to the gel //
	addFilament(tmpDNC,kAngle,visc,kT,dt,kC,fitOrder);
	// now go through the existing filaments and add crosslinks //
	Filament * fnew = filament(i);
	for(int l=0; l<i; l++) {
	  Filament * fold = filament(l);
	  VectorND sep;
	  sep = ((fold->nodes)[0])->point() - tmpDNC[0]->point();
	  _box->mapDistance(sep);
	  if(norm2(sep) <= 2.0*nodeLen*nodesPerFil) {
	    attachCrosslink(fnew,fold,kcl);
	  }
	}
      }
      _output = paraview;
      std::cout << "Set up gel with " << nFils << " filaments, " << nodesPerFil << " nodes per filament for a total of " << nFils*nodesPerFil << " nodes (self-consistency check: # of nodes in container = " << dNodes.size() << ")" << std::endl;
      std::cout << "Total # of motors = " << _motors.size() << ", total # of crosslinks = " << _crosslinks.size() << ", total # of constraints = " << _constraints.size() << "." << std::endl; 
    }

    //! virtual destructor
    virtual ~SemiflexibleGel() { 
      for(int i=0; i<_filaments.size(); i++) delete(_filaments[i]);
    }

 //    static SemiflexibleGel<N> * setupNewGel();
    
    //! Do mechanics on Filaments
    void compute( bool f0, bool f1, bool f2 );

    void printParaview(const std::string fileName) const;

    void addFilament( const DefNodeContainer & n, double kAngle, double viscosity, 
		      double kT, double dt, double kC, int fitOrder){ 
      Filament * f = new Filament( n, kAngle, viscosity, kT, dt, kC, fitOrder);
      _filaments.push_back( f );
    }

    void addFilament( Filament * f ) { _filaments.push_back( f ); }

    const FilamentContainer & filaments() const { return _filaments; }

    const Filament * filament(int a) const { return _filaments[a]; }

    Filament * filament(int a) {return _filaments[a]; }
    
    void addCrosslink( Clink * c ) { _crosslinks.push_back( c ); }

    const CrosslinkContainer & crosslinks() const { return _crosslinks; }

    const ConstraintContainer & constraints() const { return _constraints; }

    void attachCrosslink(Clink * cl) { // eventually, find nearest filaments and attempt to attach; for now, just take two filaments and attach at intersection point //
      attachCrosslink(cl, *_filaments.begin(),*(_filaments.begin()+1));
    }
  
    void attachCrosslink(Clink * cl, const Filament * f1, const Filament * f2);

    void attachCrosslink(Filament * f1, Filament * f2, double kcl);

    void addMotor(MolMot * mot) { _motors.push_back(mot); }

    void addMotor(VectorND & p, double k, double d0) {
      MolMot * mot = new Motor<N>(p,k,d0);
      _motors.push_back(mot);
    }

    void addMotor(VectorND & p, double k) {
      MolMot * mot = new Motor<N>(p,k);
      _motors.push_back(mot);
    }

    void attachMotor(MolMot * mot) { // eventually, find nearest filaments and attempt to attach; for now, just take two filaments and attach at intersection point //
      attachMotor(mot, *_filaments.begin(),*(_filaments.begin()+1));
    }

    void attachMotor(MolMot * mot, const Filament * f1, const Filament * f2);

    void attachMotor(const Filament * f1, const Filament * f2);

    const MotorContainer & motors(){ return _motors; }

    void addConstraint( Constraint * c ) { _constraints.push_back( c ); }

    void setBox( PeriodicBox * box ) { _box = box; }

    const PeriodicBox * box() const { return _box; }

    void setRNGNorm(double mn, double stddev) {
      _rngnorm = new ranlib::Normal<double>(mn,stddev);
      _rngnorm->seed((unsigned int)time(0));
    }

    void setRNGUni() {
      _rnguni = new ranlib::Uniform<double>();
      _rnguni->seed((unsigned int)time(0));
    }

  private:

    //! Filaments
    FilamentContainer 	_filaments;		

    //! Crosslinks
    CrosslinkContainer 	_crosslinks;		

    //! Constraints
    ConstraintContainer _constraints;

    //! Motors
    MotorContainer _motors;

    PeriodicBox * _box;

    CrosslinkNodeMap _crossNodeMap;

    ranlib::Normal<double> * _rngnorm;
    ranlib::Uniform<double> * _rnguni;
  };  
} // namespace voom

#include "nlSemiflexibleGel.icc"

#endif // __nlSemiflexibleGel_h__
