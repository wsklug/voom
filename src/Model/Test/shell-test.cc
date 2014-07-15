#include <string>
#include <iostream>
#include <vector>
#include <getopt.h>
#include "Node.h"
#include "SCElastic.h"
#include "FVK.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"

// #define NODENUMBER 42
// #define ELEMNUMBER 80

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{
  if( argc != 2 ) {
      cout << "Usage: shell_test modelName." << endl;
      return(0);
  }

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

  string modelName = argv[1];

  string inputFileName = modelName + ".vtk";
  // create input stream
  ifstream ifs( inputFileName.c_str());
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }
  

  // create vector of nodes
  double R = 1.0;
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;

  // find points header
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
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  Ravg /= nodes.size();
  if(verbose) cout << "Number of nodes: " <<nodes.size() << endl
		   << "Ravg = " << Ravg << endl;
  //
  // create connectivities
  while( token != "POLYGONS" ) ifs >> token;
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=0;
  ifs >> ntri;
  connectivities.reserve(ntri);
  if(verbose) cout << "Number of triangles: " <<ntri << endl;
  int ntmp=0;
  ifs >> ntmp;
  if(ntmp != 4*ntri) {
    cout << "Non-triangular polygons?" << endl;
    return 0;
  }
  for (int i = 0; i<ntri; i++){
    int tmp=0;
    ifs >> tmp;
    ifs >> tmp; c[0]=tmp;
    ifs >> tmp; c[1]=tmp;
    ifs >> tmp; c[2]=tmp;
    if(verbose) 
//       cout << setw(10) << c[0] << setw(10) << c[1] << setw(10) << c[2] 
// 	   << endl; 
    connectivities.push_back(c);
  }
  if(verbose) cout << connectivities.size() << endl;
  ifs.close();

	
  double KC = 1.0e2;
  double KG = 0.0e2;
  double C0 = 0.0;
  double E = 1.0e0;
  double nu = 0.3;
  typedef FVK matType;
  matType bilayer( KC, KG, C0, E, nu );
  std::cout << "Material has been created." << std::endl;
	
  //
  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  double penaltyVolume=1.0e3;
  double penaltyArea=1.0e2;
  double viscosity=1.0e1;
  int quadOrder=2;
  
  ifstream bdinp("body.inp");
  bdinp >> penaltyVolume >> penaltyArea >> viscosity;
  LoopShellBody<matType> bd(bilayer, connectivities, nodes, quadOrder, 
			      nBoundaries, pressure, tension, 
			      penaltyVolume, penaltyArea, viscosity,
			      LoopShell<matType>::penalty, 
			      LoopShell<matType>::penalty 
			      );
  //
  // 
  cout << "Created a body." << endl;

  // bd.printByHDS();
//   bd.compute(true,true, false);
  bd.setOutput(Body::paraview);
//   bd.printParaview("staticForce");
//   //	bd.createOpenDXData(ofn, 0);
//   //	bd.createInputFile("InputFile");
//   cout << "volume of current body = " << bd.volume() << endl;
//   //	cout << " rank test ..." << endl;
  
  //
  // create Body Container 
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

  Storage solver;
  solver.resize(model.dof());

//   model.computeAndAssemble(solver,true,true,false);

//   model.print("Init");

//   clock_t dt = 0;
//   for(int i=0; i<model.dof(); i++) {
//     clock_t t1 = clock(); 
//     for(int i=0; i<10; i++) model.computeAndAssemble(solver,true,true,false);
//     clock_t t2 = clock();
//     dt += t2-t1;
//   }
//   std::cout << 10*model.dof() << " calls of model.computeAndAssemble() took "
// 	    << dt << " clocks_t's ("
// 	    << double(dt)/CLOCKS_PER_SEC <<" seconds)."
// 	    << std::endl;

//   return 0;

  model.computeAndAssemble(solver,true,true,false);

//   std::cout << "solver.energy() = " << solver._E << std::endl
// 	    << "solver.gradient() = " << solver._DE << std::endl;

// randomly adjust positions
//  srand(time(0));
//  for(Body::NodeIterator n=bd.nodes().begin(); n!=bd.nodes().end(); n++) {
//    for(int i=0; i<(*n)->dof(); i++)
//      (*n)->addPoint(i, 0.05*rand()/RAND_MAX-0.025);
//  }
  
  bd.reduceVolume( 0.8 );

  model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
  return 0;

//   Body::ElementContainer & elements = bd.elements();
//   int eid=0;
//   blitz::Array<double,1> F(bd.nodes().size() * 3);
//   F = 0.0;
//   for(Body::ElementIterator e=elements.begin(); e!=elements.end(); e++, eid++){
//     std::cout << std::endl
// 	      << "==== element " << eid << " ====" << std::endl;
//     const ElementBase::DofIndexMap & idx = (*e)->index();
//     const ElementBase::ElementVector & f = (*e)->forces();
//     std::cout << std::endl << std::setw(8) << eid << std::setw(8) << ' ';
//     for(int a=0; a<((LoopShell<matType>*)(*e))->nodes().size(); a++) {
//       std::cout << std::setw(8) << ((LoopShell<matType>*)(*e))->nodes()[a]->id();
      //solver.gradient( idx[ei] ) += f(ei);
//       for(int ai=0; ai<3; ai++) {
// 	F(idx[3*a+ai]) += f(3*a+ai);
// 	std::cout << std::setw(12) << a 
// 		  << std::setw(12) << ((LoopShell<matType>*)(*e))->nodes()[a]->id() 
// 		  << std::setw(12) << idx[3*a+ai] 
// 		  << std::setw(12) << f(3*a+ai)  
// 		  << std::endl;
//       }
//     }
//   }
//   std::cout << F << std::endl;
//   std::cout << std::endl;

//   LoopShell<matType> * elem = 
//     (LoopShell<matType>*)((bd.elements())[0]);
//   elem->CheckConsistency();
//   return 0;

  ConjugateGradientWSK CGsolver;//(true);
  int maxIter = 1000*model.dof();
  int restartStride = model.dof()/2;
  int printStride = restartStride;
  double tol = 1.0e-12;
  double absTol = 1.0e-12;
  double tolLS = 1.0e-3;
  int maxIterLS = 20;
    
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 maxIter, restartStride, printStride, 
			 tol, absTol, tolLS, maxIterLS);
  CGsolver.setWolfeParameters(1.0e-6, 1.0e-3, 1.0e-4);
  model.print("BeforeCG");

  std::cout << "Before CG:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Area = "<< bd.area() << std::endl;
  CGsolver.solve( &model );
  model.print("AfterCG");
  std::cout << "After CG:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Area = "<< bd.area() << std::endl;
//  model.checkRank(model.dof()-6,true);

  return 0;
}


void ioSetting(int argc, char* argv[], ifstream& ifs, string& ofn)
{

  // if no arguement or too many
  if (argc < 2 || argc > 3)
    {
      cout << "Usage: ProgramName InputFileName [OutputFileName]." << endl;
      exit(0);
    }

  string inFullName = argv[1];
  string pathName;
  //	basic_string <char>::size_type sztp = inFullName.find_last_of("/");
	
  // 	if ( (sztp) != string::npos )  // no position was found
  // 		pathName = inFullName.substr(0, sztp) + "/";
  // 	else
  // 		pathName = "";

  pathName = "";
	
  // only one arguement
  if ( argc == 2 ) ofn = pathName + "output.dat";


  // two arguements
  if ( argc == 3 ) {
    ofn = argv[2];
    ofn = pathName + ofn;
  }

  // create input stream
  ifs.open( inFullName.c_str(), ios::in);
  if (!ifs)
    {
      cout << "can not open input file: " << inFullName << endl;
      exit(0);
    }
	
}
