// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug
//                University of California Los Angeles
//                   (C) 2006-2008 All Rights Reserved
//
//----------------------------------------------------------------------

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include "Node.h"
#include "CompNeoHookean.h"
#include <tvmet/Vector.h>
#include "Capsid3DBody.h"
#include "TetQuadrature.h"
#include "ShapeTet10.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "SimulatedAnnealing.h"
#include "Rigid.h" 

#if defined(_OPENMP)
#include <omp.h>
#endif

// uncomment if you would like to use for afm tip/plates
//#define AUGMENTED_LAGRANGIAN_CONTACT

using namespace tvmet;
using namespace std;
using namespace voom;

class ViscousRegularizer : public Element {

public:
  using Element::BaseNodeContainer; 
  using Element::BaseNodeIterator; 
  using Element::ConstBaseNodeIterator; 

  ViscousRegularizer( const BaseNodeContainer & nodes, double viscosity ) 
  {
    _baseNodes = nodes;
    _viscosity = viscosity;
    _energy = 0.0;
    int dof=0;
    for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      dof+=(*n)->dof();
    }
    _reference.resize(dof);
    _reference = 0.0;
    step();
  }

  void step() {
    int I=0;
    for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      for(int i=0; i<(*n)->dof(); i++,I++) {
	_reference(I) = (*n)->getPoint(i);
      }
    }
  }

  void compute(bool f0, bool f1, bool f2) {
    if(f0) _energy = 0.0;
    int I=0;
    for(BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      for(int i=0; i<(*n)->dof(); i++,I++) {
	double dx = (*n)->getPoint(i) - _reference(I);
	    if(f0) _energy += 0.5 * _viscosity * dx * dx;
	    if(f1) (*n)->addForce(i, _viscosity * dx); 
      }
    }
  }

  double viscosity() const {return _viscosity;}

private:

  double _viscosity;
  blitz::Array<double, 1> _reference;
};


int main(int argc, char* argv[])
{

#ifdef _OPENMP
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  bool verbose=true;
#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif

  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

  string objFileName = argv[1];

  ifstream ifs;
  string modelName = argv[1];
  string inputFileName = modelName + ".vtk";

  bool compressFlag=false;
  double Hfinal=0;
  bool restart=false;
  string restartFileName;
  bool boundaryConditionFlag = false;
  bool badCommandLine=false;

  for(char option; (option=getopt(argc,argv,"c:H:b:R:")) != EOF; ) 
    switch (option) {
    case 'c' :
      if( std::string(optarg) == std::string("yes") ) {
	compressFlag = true;
	std::cout << "Compress capsid." << std::endl;
      }
      break;
    case 'H' :
      Hfinal = std::atof(optarg);
      std::cout << "Hfinal: " << Hfinal << std::endl;
      break;
    case 'b' :
      if( std::string(optarg) == std::string("yes") ) {
	boundaryConditionFlag = true;
      } else if( std::string(optarg) != std::string("no") ) {
	std::cout << "Yes or no for boundary conditions.  Exiting."
		  << std::endl;
	return 0;
      }	
      std::cout << "boundaryConditionFlag: " << boundaryConditionFlag
		<< std::endl;
      break;
    case 'R' :
      restart=true;
      restartFileName = optarg;
      std::cout << "Restarting from file: " << restartFileName << std::endl;
      break;
    default :
      badCommandLine=true;
      break;
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: capsid InputFileName [-c -H finalHeight -b boundaryConditionFlag -R restartFilename]." << endl;
      return(0);
  }
  if( Hfinal <= 0.0 ) {
    std::cout << "Please specify a positive value of Hfinal." << std::endl;
    return 0;
  }

  
  // create input stream
  ifs.open( inputFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  //
  // create vector of nodes
  int dof=0;

  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  char key;  ifs>>key;
  
  // MMG, 9-18-06:  Input .vtk file containing nodes and connectivities
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;// skip number type

  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    x /= 10.0;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  cout << "Number of nodes: " <<nodes.size() << endl;

  // read in tetrahedral connectivities
  while( token != "CELLS" ) ifs >> token;
  std::vector< std::vector<int> > connectivities;
  std::vector<int> c(10);
  int ntet=0; ifs >> ntet;
  connectivities.reserve(ntet);
  cout << "Number of tetrahedrons: " << ntet << endl;
  int ntmp=0;
  ifs >> ntmp;
  if(ntmp != 11*ntet) {
    cout << "Non-quad-tetrahedral elements?" << endl;
    return 0;
  }
  for (int i = 0; i<ntet; i++){
    int tmp=0;
    ifs >> tmp;
    ifs >> tmp; c[0]=tmp;
    ifs >> tmp; c[1]=tmp;
    ifs >> tmp; c[2]=tmp;
    ifs >> tmp; c[3]=tmp;
    ifs >> tmp; c[4]=tmp;
    ifs >> tmp; c[5]=tmp;
    ifs >> tmp; c[6]=tmp;
    ifs >> tmp; c[7]=tmp;
    ifs >> tmp; c[8]=tmp;
    ifs >> tmp; c[9]=tmp;
    connectivities.push_back(c);
  }
  if(verbose) cout << connectivities.size() << endl;
  ifs.close();


  if( restart ) {
    // Load reference and current configurations from a VTK file
    std::ifstream res(restartFileName.c_str());
    std::string token;
    res >> token; 
    while( token != "POINTS" ) res >> token;
    int npts=0;
    res >> npts; 
    if(npts != nodes.size()) {
      std::cout << "Cannot restart.  " 
		<< objFileName << " and " << restartFileName 
		<< " have different numbers of nodes." << std::endl;
      return 0;
    }
    res >> token;// skip number type
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      res >> x(0) >> x(1) >> x(2);
      defNodes[i]->setPosition(x);
    }
    while( token != "displacements" ) res >> token;
    res >> token;// skip number type
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      res >> x(0) >> x(1) >> x(2);
      x += defNodes[i]->position();
      defNodes[i]->setPoint(x);
    }
        
  } 

  // create material
  double rho = 1.0;
  double E = 200.0;//1.0e3;
  double nu = 0.4;
  double viscosity=1.0e-2;

  ifstream inp("parameters.inp");
  inp >> E >> nu >> viscosity;
  if(verbose) 
    std::cout << "Input rho: " << rho << std::endl
	      << "Input E: " << E << std::endl
	      << "Input nu: " << nu << std::endl
	      << "Input viscosity: " << viscosity << std::endl;

  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  if(verbose) std::cout << "Material has been created." << std::endl;
	
  // create Body
  unsigned int quadOrder = 2;
  
  typedef Capsid3DBody<TetQuadrature,MaterialType,ShapeTet10> Capsid;
  Capsid bd(protein, connectivities, nodes, quadOrder);

  bd.setOutput(paraview);

  bd.compute(true,true,false);
  bd.printParaviewQuadTet("initial");
  std::cout << "body energy = " << bd.energy() << std::endl;

  //return 0;
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //model.checkRank(model.dof()-6,true);
  //return 0;

  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

  if( !compressFlag ) {
    std::cout << "All done.  Bye now." << std::endl;
    return 0;
  }
  
  std::cout << "Compressing capsid." << std::endl;

  double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  double Zavg = 0.0;
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
    Zmax = std::max(Zmax, Z);
    Zavg += Z;
  }
  Zavg /= defNodes.size();

  std::cout << "Zmin = " << Zmin << std::endl
	    << "Zmax = " << Zmax << std::endl
	    << "Zavg = " << Zavg << std::endl;
    

  // create indentor and plate
  double k1 = 1.0e4;
  double afmR = 10.0;
   double afmX = 0.0;
  double afmY = 0.0;
  bool frictionalContact = false;
  inp >> k1 >> afmR >> afmX >> afmY >> frictionalContact;
  if(verbose) 
    std::cout << "Input k1: " << k1 << std::endl
	      << "Input afmR: " << afmR << std::endl
	      << "Center position of afm tip (x): " << afmX << std::endl
	      << "Center position of afm tip (y): " << afmY << std::endl
	      << "Frictional Contact: " << frictionalContact << std::endl;
  
  tvmet::Vector<double,3> xc; xc = 0.0, 0.0, Zmax+afmR; 

#ifdef AUGMENTED_LAGRANGIAN_CONTACT

  RigidPlateAL * lowPlate = 0;
  
  if(verbose) std::cout << "Lower plate is a RigidPlateAL." << std::endl;
  lowPlate = new RigidPlateAL( defNodes, k1, Zmin, true);
  bd.pushBack( lowPlate );
  
#endif
  
//  RigidHemisphereContact * afm = 0;
//  if(afmR > 0.0 ) {
//    if(verbose) std::cout << "AFM is a RigidHemisphereContact." << std::endl;
//    afm = new RigidHemisphereContact( defNodes, afmR, xc );
//    bd.pushBackContact( afm );
//  }
 
 

   RigidHemisphereContact * afm = 0;
   if(afmR > 0.0 ) {
     if(verbose) std::cout << "AFM is a RigidHemisphereContact." << std::endl;
     afm = new RigidHemisphereContact( defNodes, afmR, xc, frictionalContact );
     bd.pushBackContact( afm );
   }

  // constrain top and bottom nodes

  double Rmin=std::numeric_limits<double>::max();  
  for (int a=0; a<defNodes.size(); a++) {
    double x = defNodes[a]->getPoint(0);
    double y = defNodes[a]->getPoint(1);
    double R = sqrt(x*x + y*y);
    Rmin = std::min(R,Rmin);
  }
  if(verbose) std::cout << "Rmin = " << Rmin << "." << std::endl;

  ViscousRegularizer vr(bd.nodes(), viscosity);
  bd.pushBack( &vr ); 

  int steps = 10;
  double vrTol = pgtol;
  inp >> steps;
  inp >> vrTol;
  if(verbose) 
    std::cout << "Input steps: "<< steps <<std::endl
	      << "Input vrTol: "<< vrTol <<std::endl;
  const double originalHeight = Zmax-Zmin;
  double Zbegin = Zmin;
  double Zend = Zmax-Hfinal; // /*0.5*(Zavg+Zmax);*/Zmin+20.0;
  double dZ = (Zend-Zbegin)/steps;
  double vrEnergy = vr.energy();
  double bdEnergy = bd.energy();
  
  // set up bounds for solver
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;

#ifdef AUGMENTED_LAGRANGIAN_CONTACT
  
  for(int a=0; a<defNodes.size(); a++) {
    if( afmR > 0 )      nbd(3*a+2) = 0;
    else {
      nbd(3*a+2) = 3;
      hi(3*a+2) = Zmax;
    }
  }
  
#else  

  for (int a=0; a<defNodes.size(); a++) {  
    if( afmR > 0 ) 	nbd(3*a+2) = 1;
    else {
      nbd(3*a+2) = 2;
      hi(3*a+2) = Zmax;
    }
    lo(3*a+2) = Zbegin;
  }

#endif

  solver.setBounds(nbd,lo,hi);

  if(boundaryConditionFlag) {
    // constrain apices
    for (int a=0; a<defNodes.size(); a++) {
      double x = defNodes[a]->getPoint(0);
      double y = defNodes[a]->getPoint(1);
      double R = sqrt(x*x + y*y);
      //std:: cout << R << std::endl;
      if(R < 1.0e-5) {
	std::cout << "Constraining node " << a << " at R=" << R 
		  << "." << std::endl;
	for(int i=0; i<2; i++) {
	  nbd(3*a+i) = 2;
	  hi(3*a+i) = 0.0;
	  lo(3*a+i) = 0.0;
	}
      }
    }
  }

  ofstream FvsZ;
  if(verbose) {
    FvsZ.open("FvsZ.dat");
  }

  // set xi and disp vectors for 'continuation' - a solver optimization
  std::vector< double > xin;
  std::vector< double > disp;
  
  for (int a=0; a<defNodes.size(); a++) {
    xin.push_back(defNodes[a]->getPoint(0));
    xin.push_back(defNodes[a]->getPoint(1));
    xin.push_back(defNodes[a]->getPoint(2));

    disp.push_back(0.0);
    disp.push_back(0.0);
    disp.push_back(0.0);
  }

  for(double Z = Zbegin; Z<Zend; Z+=dZ) {
    std::cout << "Z = " << Z << std::endl;

#ifdef AUGMENTED_LAGRANGIAN_CONTACT
    lowPlate->setZ( Z );
#else
    for(int a=0; a<defNodes.size(); a++) {  
      lo(3*a+2) = Z;
    }

    // move starting guess by adding prev disp to current pos, before solving.  should speed up solving.    
    for (int a=0; a<defNodes.size(); a++) {
      defNodes[a]->addPoint(0,0.99*disp[3*a]);
      defNodes[a]->addPoint(1,0.99*disp[3*a+1]);
      defNodes[a]->addPoint(2,0.99*disp[3*a+2]);      
    }
    
    if( frictionalContact ) {
      // if a node is now in contact with top or bottom, will constrain
      // both the x- and y-direction bounds before solving...but will
      // not constrain top if the AFM tip radius > 0.0
      for (int a=0; a<defNodes.size(); a++) {
	if( afmR < 0.0 ) {
	  if( defNodes[a]->getPoint(2) >= Zmax ) {	  
	    // top
	    for(int i=0; i<2; i++) {
	      nbd(3*a+i) = 2;
	      hi(3*a+i) = defNodes[a]->getPoint(i);
	      lo(3*a+i) = defNodes[a]->getPoint(i);
	    }
	  }
	} else if( defNodes[a]->getPoint(2) <= Z ) {	  
	  // bottom
	  for(int i=0; i<2; i++) {
	    nbd(3*a+i) = 2;
	    hi(3*a+i) = defNodes[a]->getPoint(i);
	    lo(3*a+i) = defNodes[a]->getPoint(i);
	  }
	}
      }
    }

    solver.setBounds(nbd,lo,hi);
#endif

    int viter=0;
    do {
      if(verbose) std::cout << "viter = " << viter << std::endl;
      viter++;
      vr.step();
      
#ifdef AUGMENTED_LAGRANGIAN_CONTACT
      double k = k1;
      double kfactor=10.0;
      int max_citer=10;
     
      for(int citer=0; citer<max_citer; citer++) {
	lowPlate->setPenaltyCoefficient(k);
	solver.solve( &model );
	if(verbose) std::cout << "citer " << citer << ": "
			      << "k = "  
			      << lowPlate->penaltyCoefficient() << ", "
			      << "penetration = " 
			      << lowPlate->penetration()
			      << std::endl; 
	double p = std::abs( lowPlate->penetration() );
	if( citer > 0 && p < 1.0e-4 ) {
	  break;
	} else {
	  k *= kfactor;
	}
      } 
#else
      solver.solve( &model );

      for (int a=0; a<defNodes.size(); a++) {
	if( defNodes[a]->getPoint(2) > Z  ) {	  
	  for(int i=0; i<2; i++) {
	    nbd(3*a+i) = 0;
	  }
	}
      }

      solver.setBounds(nbd,lo,hi);
      solver.solve( &model );

#endif
      
      vrEnergy = vr.energy();
      bdEnergy = bd.energy();

#ifdef WITH_MPI
      double temp=vrEnergy;
      MPI_Allreduce(&temp, &vrEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
      temp=bdEnergy;
      MPI_Allreduce(&temp, &bdEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
#endif      

      if(verbose)
	std::cout << "v: " << vrEnergy << std::endl
		  << "e: " << bdEnergy << std::endl;
    } while(vrEnergy > vrTol * bdEnergy);

    // calculate disp vector now that a step has solved
    // reset the xi to solved position
    for (int a=0; a<defNodes.size(); a++) {
      double xpos = defNodes[a]->getPoint(0);
      double ypos = defNodes[a]->getPoint(1);
      double zpos = defNodes[a]->getPoint(2);

      disp[3*a] = xpos - xin[3*a];
      disp[3*a+1] = ypos - xin[3*a+1];
      disp[3*a+2] = zpos - xin[3*a+2];

      xin[3*a] = defNodes[a]->getPoint(0);
      xin[3*a+1] = defNodes[a]->getPoint(1);
      xin[3*a+2] = defNodes[a]->getPoint(2);
    }

    double height = Zmax-Z;
    char name[20]; 
    sprintf(name,"disp%09.6f",originalHeight-height);
    
    // This is where Paraview file (if desired) is saved
    model.print(name);

    // add up forces on top and bottom
    double FZg = 0.0;
    double FZa = 0.0;
//     double Ztop = Zmax - 0.01*originalHeight;
//     double Zbot = Zmax - height + 0.01*originalHeight;
    double Ztop = Zmax;
    double Zbot = Zmax - height;

    for (int a=0; a<defNodes.size(); a++) {
      double F = defNodes[a]->getForce(2);
//       if(F < 0.0) 	FZa += F;
//       else if(F > 0.0) 	FZg += F;
      if( std::abs( defNodes[a]->getPoint(2) - Ztop ) < 1.0e-4 ) {
	// top
	FZa += defNodes[a]->getForce(2);
      } else if( std::abs( defNodes[a]->getPoint(2) - Zbot ) < 1.0e-4 ) {
	// bottom
	FZg += defNodes[a]->getForce(2);
      }
    }
    if(afmR > 0) FZa = afm->FZ();


#ifdef AUGMENTED_LAGRANGIAN_CONTACT
    
    FZg = lowPlate->FZ();
    
#endif

#ifdef WITH_MPI
    temp=FZg;
    MPI_Reduce(&temp, &FZg, 1, MPI_DOUBLE, 
	       MPI_SUM, 0, MPI_COMM_WORLD);
    temp=FZa;
    MPI_Reduce(&temp, &FZa, 1, MPI_DOUBLE, 
	       MPI_SUM, 0, MPI_COMM_WORLD);
#endif      
    if(verbose) {
      std::cout << "Energy = " << solver.function() << std::endl
		<< "  height = " << height
		<< "  change in height = " << originalHeight-height
		<< "  FZg = " << FZg
		<< "  FZa = " << FZa << std::endl;
      FvsZ << std::setw( 24 ) << std::setprecision(16) << originalHeight-height
	   << std::setw( 24 ) << std::setprecision(16) << FZg
	   << std::setw( 24 ) << std::setprecision(16) << FZa
	   << std::setw( 24 ) << std::setprecision(16) << solver.function()
	   << std::endl;
    }
  }
  FvsZ.close();

  std::cout << "All done.  Bye now." << std::endl;

  return 0;
}

