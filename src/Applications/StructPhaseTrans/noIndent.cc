#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "GLElastic.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"
#include "GLBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"
#include "SimulatedAnnealing.h"
#include "Lbfgsb.h"

//#define NODENUMBER 602
// #define ELEMNUMBER 80

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);


int main(int argc, char* argv[])
{

  ifstream ifs;
  string ofn;
 

  ioSetting(argc, argv, ifs, ofn);

  
  /*  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
  int dof=0;
  double eta=0.0;//reaction coordinate
  srand(time(0));    

  std::vector< NodeBase* > nodes;
  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;
    eta= (double)(rand())/RAND_MAX;
   
    XCNode<3>::Point x;
    ifs >> id >> x(0) >> x(1) >> x(2);
    //x*=6.0;
    
    NodeBase::DofIndexMap idx(4);
    
   for(int j=0; j<4; j++) idx[j]=dof++;
   nodes.push_back(new XCNode<3>(id,idx,x,eta));    

}
    


  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; i < ELEMNUMBER; i++){
    ifs >> c[0];//c[0]-=1;
    ifs >> c[1];//c[2]-=1;
    ifs >> c[2];//c[1]-=1;

    connectivities.push_back(c);
  }

  ifs.close();
  */

  std::vector< NodeBase* > nodes;
  std::vector< XCNode<3>* > xcNodes;
  double Ravg = 0.0;
  int dof = 0;

  // find points header
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  nodes.reserve(npts);
  xcNodes.reserve(npts);
  ifs >> token;// skip number type

  double eta = 0.0; //reaction coordinate
  srand(time(0));  
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    XCNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    double r=tvmet::norm2(x);
    //x*=rho;
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(4);

    for(int j=0; j<4; j++) idx[j]=dof++;
    //random distribution of eta
//     if (r>=0.95)
//       eta=0.0;
//     else
//       eta=1.0;
    eta= (double)(rand())/RAND_MAX;
    XCNode<3>* n = new XCNode<3>(id,idx,x,eta);
    nodes.push_back( n );
    xcNodes.push_back(n);
  }
  Ravg /= nodes.size();
  cout << "Number of nodes: " <<nodes.size() << endl
       << "Ravg = " << Ravg << endl;

  // read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=0; ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "Number of triangles: " <<ntri << endl;
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

    connectivities.push_back(c);
  }

  ifs.close();        

  double kC = 1.0e0;
  double kG = 0.0e0;
  double C0 = 0.0;
  double rho = 6.0;    //rho=R/delta, delta = 2 nm is size of a protein
  double gamma = 100.0;//FvK number
  double Y  = 1.0;    
  double nu = 1.0/3.0;
  double DeltaEa = 0.01;
  double epsilon = 0.01;
  double DeltaEcPrime = 0.01; //DeltaEc'
  double Gamma = 1.0;
  double viscosity = 0.0;

  
  ifstream bdinp("constants.inp");
  bdinp >> rho >> gamma  >> nu >> epsilon >> Gamma >> DeltaEcPrime >> DeltaEa;
  std::cout << "rho ="     << rho << std::endl
            << "gamma ="   << gamma  <<  std::endl
            << "nu ="      << nu <<  std::endl
	    << "epsilon =" << epsilon    << std::endl
	    << "Gamma ="     << Gamma        << std::endl
	    << "DeltaEc' ="<< DeltaEcPrime   << std::endl
	    << "DeltaEa ="  << DeltaEa     << std::endl;

  double bb = -2.0-2.0*DeltaEcPrime;
  double cc = 1.0+3.0*DeltaEcPrime;
  double etaPrime = 0.5+1.5*DeltaEcPrime;
  double psiZero = DeltaEa/etaPrime/etaPrime/(etaPrime*etaPrime + bb*etaPrime + cc);
  double DeltaEc = DeltaEcPrime * psiZero;

  std::cout << "DeltaEc ="<< DeltaEc <<std::endl
	    << "psiZero ="<< psiZero <<std::endl;
  Y = gamma/rho/rho;
  std::cout << "Y ="<<Y<<std::endl;

  GLElastic bilayer(kC, kG, C0, Y, nu,  DeltaEa, epsilon, Gamma, DeltaEcPrime);
  //SCElastic bilayer(KC, KG, C0);
  std::cout << "GL Material has been created." << std::endl;
  //std::cout << "SCElastic Material has been created." << std::endl;	
  //	

  //rescale the size
  for(int i=0; i<xcNodes.size(); i++) {
    XCNode<3>::Point x;
    x = xcNodes[i]->point();
    x *= rho/Ravg;
    xcNodes[i]->setPoint(x);
    xcNodes[i]->setPosition(x);
  }

  //Start from a sphere
//   for (vector<NodeBase*>::iterator np = nodes.begin(); np != nodes.end(); np ++){
//     XCNode<3>::Point x;
//     for(int i=0; i<3;i++)
//       x(i)=(*np)->getPoint(i);
//     double Rfactor=tvmet::norm2(x);
//     x*=rho/Rfactor;
//     for(int i=0; i<3;i++)
//       (*np)->setPoint(i,x(i));
//   }

  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  int quadOrder=2;
  double penaltyVolume=0.0;
  double penaltyArea=0.0;        
   GLBody<GLElastic> bd(bilayer, connectivities, nodes, quadOrder, nBoundaries, 
                        pressure, tension, penaltyVolume, penaltyArea, viscosity); 

  //LoopShellBody<SCElastic> bd(bilayer, connectivities, nodes, ngp, 
  //		      nBoundaries, pressure, penaltyVolume, penaltyArea);
  cout << "Created a body." << endl
       << "  # nodes: " << bd.nodes().size()<<std::endl
       << "  # shell elements: " << bd.elements().size() << std::endl;
  
  bd.setOutput(Body::paraview);
  bd.compute(false,false,false);


  std::cout << "  # all elements: " << bd.elements().size() << std::endl;

  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

  //  model.checkConsistency(true,false);
  //return 0;

  int m=5;
  double factr=0.1;//1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  double pentol=1.0e-4;
  //bool boundaryConditionFlag = true;

  //ifstream lbfgsbinp("lbfgsb.inp");
  //lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  //if(verbose) 
   std::cout  << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);
  
//   SimulatedAnnealing solver;

//   unsigned nSteps=1000; 

//   double finalTratio=1.0e-6;
//   unsigned printStride=5;
//   double T01=1.0;
//   double T02=0.01;

//   solver.setParameters(voom::SimulatedAnnealing::EXPONENTIAL,
// 		       nSteps, T01, T02, finalTratio, printStride);

  
  model.print("BeforeBFGS");
  
  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);

  std::cout << "Before BFGS:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Reduced Volume = " << vred << std::endl
	    << "Energy = " << solver.function() << std::endl
	    << "strainEnergy = " << bd.totalStrainEnergy() << std::endl;

  //model.checkConsistency(true,false);
  //return 0;
  
  solver.solve( &model );      

  //int g=gamma;
  char name[20];sprintf(name, "CapsidPsi6-%f", Gamma);
  model.print(name);
  
  v = bd.volume();
  a = bd.area(); 
  vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);

    std::cout  << "Volume = "<< bd.volume() << std::endl
	       << "Area = "<< bd.area() << std::endl
	       << "Reduced Volume = " << vred << std::endl
	       << "Energy = " << solver.function() << std::endl
	       << "strainEnergy = " << bd.totalStrainEnergy() << std::endl;



    Ravg=0.0;

  for (vector<NodeBase*>::iterator np = nodes.begin(); np != nodes.end(); np ++){
    XCNode<3>::Point x;
    for(int i=0; i<3;i++)
      x(i)=(*np)->getPoint(i);
    Ravg+=tvmet::norm2(x);
  } 

  Ravg/=nodes.size();
  std::cout  << "Ravg = "<< Ravg << std::endl;

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
