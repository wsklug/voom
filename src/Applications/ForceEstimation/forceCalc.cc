#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "SCElastic.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"
#include "SimulatedAnnealing.h"

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);


int main(int argc, char* argv[])
{
  ifstream ifs;
  string ofn;
 
  ioSetting(argc, argv, ifs, ofn);
	
  //
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  ofstream inputVTK("input.vtk");
  inputVTK <<"# vtk DataFile Version 2.0\nTest example\nASCII\nDATASET POLYDATA"
	   << std::endl;
  char key;
  ifs>>key;
  for ( int i = 0; key=='v'; i++, ifs>>key){
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
//     cout << setw(12) << id 
// 	 << setw(20) << x(0)
// 	 << setw(20) << x(1)
// 	 << setw(20) << x(2) << endl;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    nodes.push_back(new DeformationNode<3>(id,idx,x));
  }
  double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  cout << nodes.size() << endl;
  inputVTK << "POINTS " << nodes.size() << " double" << std::endl;

  for(int i=0; i<nodes.size(); i++) {
//     cout << setw(12) << nodes[i]->id()
    inputVTK << setw(20) << nodes[i]->getPoint(0)
	     << setw(20) << nodes[i]->getPoint(1)
	     << setw(20) << nodes[i]->getPoint(2) 
	     << endl;
    Zmin = std::min(Zmin,nodes[i]->getPoint(2));
    Zmax = std::max(Zmax,nodes[i]->getPoint(2));
  }
  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; key=='f'; i++){
    int tmp=0;
    ifs >> tmp; c[0]=tmp-1;
    ifs >> tmp; c[1]=tmp-1;
    ifs >> tmp; c[2]=tmp-1;
    connectivities.push_back(c);
//     cout << setw(12) << connectivities[i][0] 
// 	 << setw(12) << connectivities[i][1] 
// 	 << setw(12) << connectivities[i][2]
// 	 << endl;    
    if(!(ifs>>key)) break;
  }
  cout << connectivities.size() << endl;
  inputVTK << "POLYGONS  " << connectivities.size() 
	   << "  " << 4*connectivities.size() << endl;
  for(int i=0; i<connectivities.size(); i++) {
//     cout << setw(12) << connectivities[i][0] 
    inputVTK << 3 
	     << setw(12) << connectivities[i][0]
	     << setw(12) << connectivities[i][1] 
	     << setw(12) << connectivities[i][2]
	     << endl;
  }
  inputVTK.close();
  ifs.close();
	
  double KC = 1.0e0;
  double KG = 0.0;
  double C0 = 0.0;
  SCElastic bilayer( KC, KG, C0 );
  std::cout << "SCElastic Material has been created." << std::endl;
	
  //
  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  double penaltyVolume=1.0e7;
  double penaltyArea=1.0e6;
  double viscosity=1.0e-2;
  int quadOrder = 2;
  
  ifstream inp("parameters.inp");
  inp >> pressure >> tension;
  std::cout << "Input pressure: " << pressure << std::endl
	    << "Input tension:  " << tension << std::endl;
  typedef LoopShellBody<SCElastic> LSB;
  LSB bd(bilayer, connectivities, nodes, quadOrder, 
	 nBoundaries, pressure, tension,
	 penaltyVolume, penaltyArea, viscosity,
	 LoopShell<SCElastic>::multiplier, 
	 LoopShell<SCElastic>::multiplier 
	 );
  cout << "Created a body." << endl
       << " volume = " << bd.volume() << endl;

  bd.setOutput(Body::paraview);

  bd.compute(false,true,false);
  ofstream FvsZ("FvsZ");
  std::cout << std::setw(18) << "z-height" << std::setw(18) << "Bead Force" 
	    << std::endl;
  double Zlen = Zmax-Zmin;
  double Zstep = Zlen/100.0;
  for(double slice=Zmin; slice<=Zmax; slice+=Zstep) {
    double beadForce=0.0;
    for(LSB::ConstNodeIterator n=bd.nodes().begin(); n!=bd.nodes().end(); n++) {
      LSB::FeNode_t * sn = dynamic_cast<LSB::FeNode_t*>(*n);
      if( sn ) {
	if( sn->getPoint(2) >= slice ) beadForce += sn->getForce(2);
      }
    }
    std::cout << std::setw(18) << slice 
	      << std::setw(18) << beadForce << std::endl;
    FvsZ << std::setw(18) << slice 
	 << std::setw(18) << beadForce << std::endl;
  }

  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

  string fname = argv[1];
  fname += ".shape";
  model.print(fname);

  return 0;

//   model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
	
  ConjugateGradientWSK CGsolver;//(true);
  int maxIter = 1000*model.dof();
  int restartStride = 2*model.dof();
  int printStride = 10*restartStride;
  double tol = 1.0e-4;//8;
  double absTol = 1.0e-4;//8;
  double tolLS = 1.0e-6;
  int maxIterLS = 20;
    
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 maxIter, restartStride, printStride, 
			 tol, absTol, tolLS, maxIterLS);
  CGsolver.setWolfeParameters(1.0e-6, 1.0e-3, 1.0e-4);

  model.print("BeforeCG");
  
  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = bd.constraintVolume(); 
  double A = bd.constraintArea();
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
  std::cout << "Before CG:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Ref. Volume = "<< bd.constraintVolume() << std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Ref. Area = "<< bd.constraintArea() << std::endl
	    << "Reduced Volume = " << vred << std::endl
	    << "Ref. Reduced Volume = " << Vred << std::endl
	    << "Energy = " << CGsolver.function() << std::endl;
//   SimulatedAnnealing SAsolver;
//   const unsigned nSteps=100*model.dof(); 
//   double T01=1.0e3;
//   double T02=1.0e-1;
//   cout << "T01: ";
//   cin >> T01;
//   cout << "T02: ";
//   cin >> T02;
//   printStride=10*model.dof();
//   SAsolver.setParameters( SimulatedAnnealing::EXPONENTIAL,
// 			  nSteps, T01, T02, printStride);
//   SAsolver.solve( &model );
//   v = bd.volume();
//   a = bd.area(); 
//   vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
//   V = bd.constraintVolume();
//   A = bd.constraintArea(); 
//   Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
//   std::cout << "Volume = "<< bd.volume() << std::endl
// 	    << "Ref. Volume = "<< bd.constraintVolume() << std::endl
// 	    << "Area = "<< bd.area() << std::endl
// 	    << "Ref. Area = "<< bd.constraintArea() << std::endl
// 	    << "Reduced Volume = " << vred << std::endl
// 	    << "Ref. Reduced Volume = " << Vred << std::endl
// 	    << "Energy = " << CGsolver.function() << std::endl;
//   return 0;
  for(double nu=1.0; nu > 0.5; nu-=0.1) {
    for(int iter=0; iter<5; iter++) {
      bd.reduceVolume(nu);
      bd.resetReference();
      CGsolver.solve( &model );
    }
    char name[20]; sprintf(name,"nu%f",nu);
    model.print(name);
    v = bd.volume();
    a = bd.area(); 
    vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
    V = bd.constraintVolume();
    A = bd.constraintArea(); 
    Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
    std::cout << "nu = "<< nu << std::endl
	      << "Volume = "<< bd.volume() << std::endl
	      << "Ref. Volume = "<< bd.constraintVolume() << std::endl
	      << "Area = "<< bd.area() << std::endl
	      << "Ref. Area = "<< bd.constraintArea() << std::endl
	      << "Reduced Volume = " << vred << std::endl
	      << "Ref. Reduced Volume = " << Vred << std::endl
	      << "Energy = " << CGsolver.function() << std::endl;
  }
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
