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

#if !defined(__SemiflexibleGel_h__)
#define __SemiflexibleGel_h__

#include<blitz/array.h>
#include<random/exponential.h>
#include<vector>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <string>
#include <fstream>
#include "Body.h"

#include "Spring.h"
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
#include "Grid.h"
#include "TwoBodyPotential.h"
#include "IntersectionFinder.h"
#include "PeriodicTie.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "PinchForce.h"
#include "NematicProbTable.h"
#include "AffinityMeasure.h"

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

    typedef DeformationNode<N> BaseDefNode;

    typedef Spring<N> Bond;
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
	
    //Mo	
    typedef Crosslink<N> Clink;
    typedef std::vector< Clink* > CrosslinkContainer;
    typedef typename CrosslinkContainer::iterator CrosslinkIterator;
    typedef typename CrosslinkContainer::const_iterator ConstCrosslinkIterator;
  
    typedef Motor<N> MolMot;
    typedef std::vector< MolMot* > MotorContainer;
    typedef typename MotorContainer::iterator MotorIterator;
    typedef typename MotorContainer::const_iterator ConstMotorIterator;

    typedef PinchForce<N> Pinch;
    typedef std::vector< Pinch* > PinchContainer;
    typedef typename PinchContainer::iterator PinchIterator;
    typedef typename PinchContainer::const_iterator ConstPinchIterator;

    typedef typename std::set<DefNode*> PinchNodeSet;
    typedef typename PinchNodeSet::iterator PinchNodeIterator;
		
    typedef tvmet::Vector<double,N> VectorND;
    typedef tvmet::Matrix<double,N,N> TensorND;

    typedef typename std::map< DefNode*, DefNode* > CrosslinkNodeMap;
    typedef typename CrosslinkNodeMap::iterator CLNMiter;

    typedef typename std::map< double, int > CrosslinkDistFreq;
    typedef typename CrosslinkDistFreq::iterator ClDFIter;
    typedef std::pair< int, int > intPair;
    typedef typename std::multimap< DefNode*, int >::iterator SNiter;

    typedef std::pair< double, double > doublePair;
    typedef typename std::vector< doublePair > doublePairContainer;

    typedef std::pair< double, doublePair > doublePairWErrors;
    typedef typename std::vector< doublePairWErrors > doublePairWErrorsContainer;

    typedef typename std::map< std::string, std::string > PropertyList;
    typedef typename PropertyList::iterator PropertyIterator;
    typedef typename PropertyList::const_iterator ConstPropertyIterator;

    typedef typename std::vector< std::pair<VectorND,TensorND> > StrainField;
    typedef typename StrainField::iterator StrainFieldIterator;
	
    struct Filament {
      //! Construct with linear springs
      Filament( const DefNodeContainer & n, double kBond, double kAngle, 
		double viscosity, double kT, double dt);

      //! Construct with nonlinear (entropic) springs
      Filament( const DefNodeContainer & n, double kAngle, double viscosity, 
		double kT, double dt, double kC, int fitOrder, double maxForce);

      Filament(const DefNodeContainer & n, double kappa, double mu, double viscosity, double kT, double dt, double minLength);

      ~Filament();

      const VectorND & point();

      DefNodeContainer nodes;
      BondContainer    bonds;
      AngleContainer   angles;
      RodContainer     rods;
      std::vector< double > clinks;

      VectorND pt;
    };

    typedef std::vector< Filament* > FilamentContainer;
    typedef typename FilamentContainer::iterator FilamentIterator;
    typedef typename FilamentContainer::const_iterator ConstFilamentIterator;
    
    typedef Grid<Filament,Filament,N> FilGrid;
    typedef Grid<DefNode,BaseDefNode,N> NodeGrid;
    typedef Grid<AffinityElement,AffinityElement,N> AffElementGrid;

    struct TempCrosslink {
      //VectorND location;
      int baseFil;
      std::map< int, VectorND > otherFils;
      bool active;
      DefNode* clNode;
      PeriodicTie<N> * ptie;
    };
 
    struct TempFilament {
      VectorND start;
      VectorND end;
      std::vector< std::pair<VectorND,VectorND> > filSegs;
      std::multimap<double,TempCrosslink *> crossFils;
    };

    struct TempBox {
      std::vector<TempCrosslink *> boxCLs;
    };

    typedef typename std::vector<TempFilament *> TempFilamentContainer;

    //! Default Constructor
    SemiflexibleGel() {
      _output=paraview;
    }

    //! Another Constructor
    SemiflexibleGel( DefNodeContainer & dNodes, PeriodicBox * box, 
		     double filDens, int nodesPerFil, double nodeLen, 
		     const string & bondType, bool cutOffEnds, const PropertyList & properties );    

    //! Constructor that reads in data from a file
    SemiflexibleGel(std::string fileName, DefNodeContainer & nodes, std::string bondType, bool cutOffEnds, const PropertyList & properties, double minLength);

    //! Constructor that implements adaptive meshing
    SemiflexibleGel(DefNodeContainer & dNodes, PeriodicBox * box, double filDens, double filLength, const string & bondType, bool cutOffEnds, double minLength, const PropertyList & properties);

    //! Constructor that reads in data from a file and adaptively meshes
    SemiflexibleGel(std::string fileName, DefNodeContainer & dNodes, std::string bondType, bool cutOffEnds, double minLength, const PropertyList & properties);

    //! virtual destructor
    virtual ~SemiflexibleGel() { 
      for(int i=0; i<_filaments.size(); i++) delete(_filaments[i]);
    }

    //! reset nodal positions to initial positions, reset energy to 0
    void resetGel() {
      for(FilamentIterator f=_filaments.begin(); f!=_filaments.end(); f++) {
	for(DefNodeIterator n=f->nodes.begin(); n!=f->nodes.end(); n++) {
	  (*n)->setPoint((*n)->position());
	}
      }
      _energy = 0.0;
    }
    
    void removePrestress();
    
    int removeCrosslinks(TempFilamentContainer& tmpFils, double minLength);

    int removeLongFilCrosslinks(TempFilamentContainer& tmpFils, double L);

    double computeLongFilCrosslinks(TempFilamentContainer & tmpFils, double L);    

    int collapseCrosslinks(TempFilamentContainer & tmpFils, double minLength);

    //! Do mechanics on Filaments
    void compute( bool f0, bool f1, bool f2 );

    void printParaview(const std::string fileName) const;
    
    void storeGel(std::string fileName);

    void storeSparseGel(std::string fileName, TempFilamentContainer & tmpFils);

    void addFilament( const DefNodeContainer & n, double kBond, double kAngle,
		      double viscosity, double kT, double dt ) { 
      Filament * f = new Filament( n, kBond, kAngle, viscosity, kT, dt );
      _filaments.push_back( f );
    }

    //! Add a filament with nonlinear (entropic) springs
    void addFilament( const DefNodeContainer & n, double kAngle, double viscosity, 
		      double kT, double dt, double kC, int fitOrder, double maxForce){ 
      Filament * f = new Filament( n, kAngle, viscosity, kT, dt, kC, fitOrder, maxForce);
      _filaments.push_back( f );
    }

    void addFilament(const DefNodeContainer & n, double kappa, double mu, double viscosity, double kT, double dt, double minLength) {
      Filament * f = new Filament(n,kappa,mu,viscosity,kT,dt,minLength);
      _filaments.push_back(f);
    }

    void addFilament( Filament * f ) { _filaments.push_back( f ); }

    const FilamentContainer & filaments() const { return _filaments; }

    const Filament * filament(int a) const { return _filaments[a]; }

    Filament * filament(int a) {return _filaments[a]; }
    
    void moveCLNodes(Filament * f);

    void setGrid(FilGrid * g) { _grid = g; }

    void addTwoBodyPotential(TwoBodyPotential * tbp) { _tbp.push_back(tbp); }

    void addCrosslink( Clink * c ) { _crosslinks.push_back( c ); }

    const CrosslinkContainer & crosslinks() const { return _crosslinks; }

    const ConstraintContainer & constraints() const { return _constraints; }

    void attachCrosslink(Clink * cl) { // eventually, find nearest filaments and attempt to attach; for now, just take two filaments and attach at intersection point //
      attachCrosslink(cl, *_filaments.begin(),*(_filaments.begin()+1));
    }
  
    void attachCrosslink(Clink * cl, const Filament * f1, const Filament * f2);

    bool attachCrosslink(Filament * f1, Filament * f2, double kcl, double relax);

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

    void addPinches(double pinchDensity, double a, double tol, double f0);
    
    void addPinches(int nPinches, double a, double tol, double f0);

    void addPinch(DefNode * n1, DefNode * n2, double f0);

    void addPinch(double a, double tol, double f0);

    //void addPinch(double f0, bool springy, DefNodeContainer & dNodes, double kBond, double kAngle, double visc, double kT, double dt, double kcl);

    void setBox( PeriodicBox * box ) { _box = box; }

    PeriodicBox * box() { return _box; }

    double getMeanCLsep() { return _meanCLsep; }
    
    double getNematicOP() { return _nematicOP; }

    VectorND & getNemDirector() { return _nemDirector; }

    double getMeanFilLen() { return _meanFilLen; }

    double getMeanFilStreStiff() { return _meanFilStreStiff; }

    double getMeanFilBendStiff() { return _meanFilBendStiff; }

    double crosslinkenergy();

    double bendingenergy();

    double stretchingenergy();

    double parallelenergy();
    
    double perpenergy();

    CrosslinkDistFreq & getCrossDistro() {
      CrosslinkDistFreq & cdf = _crossDistFreqs;
      return cdf;
    }

    std::map< double, int > & getLengthDistro() {
      std::map< double, int > & ld = _filLenFreqs;
      return ld;
    }

    std::map< double, int > & getStreStiffDistro() {
      std::map< double, int > & ssd = _filStreStiffFreqs;
      return ssd;
    }

    std::map< double, int > & getBendStiffDistro() {
      std::map< double, int > & bsd = _filBendStiffFreqs;
      return bsd;
    }

    std::map< double, int > & getNematicDistro() {
      std::map< double, int > & nd = _nematicFreqs;
      return nd;
    }

    std::map< doublePair, doublePair > getAngularEnergyDistro();

    std::multimap< double, std::vector<double> > getDensityEnergyDistro(double scale);

    void cutOffEndsandCCD(double kcl, DefNodeContainer & dNodes);
    
    void computeCrossDistro(double kcl);

    void computeNematicDistro(double nemAngle);

    void computeFilLenDistro();

    void computeFilStreStiffDistro();

    void computeFilBendStiffDistro();

    void computeNematicOP();

    std::vector< std::pair<double,double> > computeNemCorrelations(double minSep, double step, double tol);

    void printAngles(std::string & angleFile);	

    bool isSlave(DefNode* node) {
      if(_crossNodeMap.find(node) == _crossNodeMap.end()) return false;
      else if(_crossNodeMap[node] == node) return false;
      else return true;
    }

    bool isMaster(DefNode* node) {
      if(_crossNodeMap.find(node) == _crossNodeMap.end()) return false;
      else if(_crossNodeMap[node] != node) return false;
      else return true;
    }

    bool areLinked(DefNode* node1, DefNode* node2) {
      if(_crossNodeMap.find(node1) == _crossNodeMap.end() || _crossNodeMap.find(node2) == _crossNodeMap.end()) return false;
      else {
	DefNode* n1mast = _crossNodeMap.find(node1)->second;
	DefNode* n2mast = _crossNodeMap.find(node2)->second;
	if(node1 == n2mast || node2 == n1mast) return true;
	else return false;
      }
    }
    
    bool checkCrosslinks();

    void checkParallelForces();

    void printCLFilDist();
    
    void printInitialBends();

    void printBigBends();

    doublePairContainer affineMeasurement(double minLength, double stepSize, double maxLength, double shear, std::string measureType);

    doublePairWErrorsContainer affineMeasurementHeadLevine(double minLength, double stepSize, double maxLength, double shear, double smallestFil, double largestFil, bool getAngularDist);

    doublePairWErrorsContainer affineMeasurementHeadLevineInterpolated(double minLength, double stepSize, double maxLength, double shear, double largestFil);

    void affineBoxesMeasurement(VectorND & boxsize, double pairDist, double shear, double largestFil, std::string fileName, bool doCorrTest);

    doublePairContainer energyCorrelationFunction(double boxs, double maxlen);

    void computeNonaffinityLengthDensityCorrelation(double boxsize, double maxdist, double maxFL, double shear);

    doublePairContainer computeCorrelationFunction(std::map<DefNode *, double> & dataPts, double minSep, double maxSep, double step, double tol);

    doublePairContainer computeCrossCorrelationFunction(std::map< DefNode *, std::pair<double,double> > & dataPts, double minSep, double maxSep, double step, double tol);

    double computeCorrelation(std::vector<double> & dat1, std::vector<double> & dat2);

    void computeCrossCorrelations(double len, double shear, std::string & fileName);

    double computeBucklingEnergy(double shear, double kap, double mu);

    void printLongShortStats(double shear, double L, double longshortratio);

  private:

    //! Filaments
    FilamentContainer 	_filaments;		

    //! Crosslinks
    CrosslinkContainer 	_crosslinks;		

    //! Constraints
    ConstraintContainer _constraints;

    //! Motors
    MotorContainer _motors;

    PinchContainer _pinches;

    PeriodicBox * _box;

    CrosslinkNodeMap _crossNodeMap;

    PinchNodeSet _pinchNodes;

    std::map< DefNode*, int > _nSlavesMap;

    CrosslinkDistFreq _crossDistFreqs;

    std::set<DefNode*> _crosslinkNodes;

    std::map< double, int > _filLenFreqs;

    std::map< double, int > _filStreStiffFreqs;

    std::map< double, int > _filBendStiffFreqs;

    std::map< double, int > _nematicFreqs;

    double _meanCLsep;
    
    double _nematicOP;

    VectorND _nemDirector;

    double _meanFilLen;

    double _meanFilStreStiff;

    double _meanFilBendStiff;

    std::vector<TwoBodyPotential*> _tbp;
    
    FilGrid * _grid;

  };  
} // namespace voom

#include "SemiflexibleGel.icc"

#endif // __SemiflexibleGel_h__
