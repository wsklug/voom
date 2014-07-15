#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
//#include "TwoPhaseElastic.h"
#include "EvansElastic.h"
#include <tvmet/Vector.h>
#include <fstream>
//#include "TwoPhaseBody.h"
#include "LoopShellBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"
#include "SimulatedAnnealing.h"
#include "Lbfgsb.h"

// #define NODENUMBER 42
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

	
  /*  
  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
   srand(time(0));
   int dof=0;
   double con=0.0; //concentration

  std::vector< NodeBase* > nodes;
  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;
    XCNode<3>::Point x;
    con = (double)(rand())/RAND_MAX;//random concentrations

    ifs >> id >> x(0) >> x(1) >> x(2);
    x(2)*=1.2;
    cout << setw(12) << id 
	 << setw(20) << x(0)
	 << setw(20) << x(1)
	 << setw(20) << x(2) 
	 << setw(20) << con << endl;

    NodeBase::DofIndexMap idx(4);
    
   for(int j=0; j<4; j++) idx[j]=dof++;
   nodes.push_back(new XCNode<3>(id,idx,x, con));    

}
    

  cout << nodes.size() << endl;
  for(int i=0; i<NODENUMBER; i++) {

    cout << setw(12) << nodes[i]->id()
	 << setw(20) << nodes[i]->getPoint(0)
	 << setw(20) << nodes[i]->getPoint(1)
	 << setw(20) << nodes[i]->getPoint(2) 
	 << setw(20) << nodes[i]->getPoint(3)
	 << endl;
  }

  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; i < ELEMNUMBER; i++){
    ifs >> c[0];
    ifs >> c[1];
    ifs >> c[2];
    connectivities.push_back(c);
  }
  cout << connectivities.size() << endl;
  for(int i=0; i<ELEMNUMBER; i++) {
    cout << setw(12) << connectivities[i][0] 
	 << setw(12) << connectivities[i][1] 
	 << setw(12) << connectivities[i][2]
	 << endl;
  }
  ifs.close();
        
  
  double kC1 = 2.0e0;
  double kG1 = 2.0e0;
  double C01 = 0.0;
  double kC2 = 2.0e0;
  double kG2 = 2.0e0;
  double C02 = 0.0;
  double mu = 3.5e0;
  double kS = 2.0e0;
  double epsilon = 0.01;
  double h = 1.0;
  double viscosity = 0.0e0;
  

  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  double chemicalTension=0.0;
  int quadOrder=1;
  double penaltyVolume=1.0e4;
  double penaltyArea=1.0e4;
  double penaltyAreaOne=1.0e4;
  
  ifstream bdinp("constants1.inp");
  bdinp >> mu >> kS >>penaltyVolume >> penaltyArea >> penaltyAreaOne >> epsilon >> h >> viscosity;
  std::cout << "mu ="  << mu <<  std::endl
            << "kS ="  << kS <<  std::endl
            << "pV ="  << penaltyVolume << std::endl
            << "pA ="  << penaltyArea   << std::endl
	    << "pA1="  << penaltyAreaOne<< std::endl
            << "epsilon = " << epsilon  << std::endl
	    << "h ="   << h             << std::endl
	    << "viscosity=" << viscosity<< std::endl;
       
  
  TwoPhaseElastic bilayer(kC1, kG1, C01, kC2, kG2, C02, mu, kS, epsilon, h);         
  std::cout << "TwoPhaseElastic Material has been created." << std::endl;

  TwoPhaseBody<TwoPhaseElastic> bd(bilayer, connectivities, nodes, quadOrder, 
				   nBoundaries,  pressure, tension, chemicalTension, 
				   penaltyVolume, penaltyArea, penaltyAreaOne, viscosity);

  cout << "Created a body." << endl
       << "  # nodes: " << bd.nodes().size()<<std::endl
       << "  # shell elements: " << bd.elements().size() << std::endl;
  
  bd.setOutput(Body::paraview);
  bd.compute(false,false,false);
	    cout<< "Volume = "<< bd.volume() << std::endl	        
	        << "Area = "<< bd.area() << std::endl;
  
  */
  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
  int dof=0;

  std::vector< NodeBase* > nodes;
  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;
    DeformationNode<3>::Point x;
    ifs >> id >> x(0) >> x(1) >> x(2);
    //x(2)*=0.8;
    cout << setw(12) << id 
	 << setw(20) << x(0)
	 << setw(20) << x(1)
	 << setw(20) << x(2) << endl;
    NodeBase::DofIndexMap idx(3);
    
   for(int j=0; j<3; j++) idx[j]=dof++;
   nodes.push_back(new DeformationNode<3>(id,idx,x));    

}
  
  double beta =1.3;//1.25;//1.25;
  for(int i=0; i<NODENUMBER; i++) {
    double zz = beta * nodes[i]->getPoint(2);
    nodes[i]->setPoint(2,zz);   
  }

  cout << nodes.size() << endl;
  for(int i=0; i<NODENUMBER; i++) {
    for(int j=0; j<3; j++) 
    cout << setw(12) << nodes[i]->id()
	 << setw(20) << nodes[i]->getPoint(0)
	 << setw(20) << nodes[i]->getPoint(1)
	 << setw(20) << nodes[i]->getPoint(2) 
	 << endl;
  }

  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; i < ELEMNUMBER; i++){
    ifs >> c[0];
    ifs >> c[1];
    ifs >> c[2];
    connectivities.push_back(c);
  }
  cout << connectivities.size() << endl;
  for(int i=0; i<ELEMNUMBER; i++) {
    cout << setw(12) << connectivities[i][0] 
	 << setw(12) << connectivities[i][1] 
	 << setw(12) << connectivities[i][2]
	 << endl;
  }
  ifs.close();

  double KC = 2.0e0;
  double KG = 0.0e0;
  double C0 = 0.0;
  double mu = 3.5e0;
  double KS = 2.0e0;
  double viscosity = 0.0e0;
  

  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  int quadOrder=2;
  double penaltyVolume=1.0e4;
  double penaltyArea=1.0e4;
  
  ifstream bdinp("constants.inp");
  bdinp >> mu >> KS >>penaltyVolume >> penaltyArea >>viscosity;
  std::cout << "mu ="  << mu <<  std::endl
            << "KS ="  << KS <<  std::endl
            << "pV ="  << penaltyVolume << std::endl
            << "pA ="  << penaltyArea   << std::endl
            << "viscosity =" << viscosity  << std::endl;
  
          
  EvansElastic bilayer( KC, KG, C0, mu, KS );
  //SCElastic bilayer(KC, KG, C0);
  std::cout << "EvansElastic Material has been created." << std::endl;
  //std::cout << "SCElastic Material has been created." << std::endl;	
  //

  LoopShellBody<EvansElastic> bd(bilayer, connectivities, nodes, quadOrder, 
				 nBoundaries, pressure, tension, penaltyVolume, penaltyArea, viscosity); 


  //LoopShellBody<SCElastic> bd(bilayer, connectivities, nodes, ngp, 
  //		      nBoundaries, pressure, penaltyVolume, penaltyArea);
  cout << "Created a body." << endl
       << "  # nodes: " << bd.nodes().size()<<std::endl
       << "  # shell elements: " << bd.elements().size() << std::endl;


  bd.setOutput(Body::paraview);
  bd.compute(false,false,false);

  srand(time(0));
   for(int i=0; i<NODENUMBER; i++) {
     for(int j=0; j<3; j++) 
       nodes[i]->addPoint(j,0.01*((double)(rand())/RAND_MAX-0.5) );
   }
  

  cout<< "Volume = "<< bd.volume() << std::endl	        
      << "Area = "<< bd.area() << std::endl;	
  
  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = -1;//0;
  double pentol=1.0e-4;
  //ifstream lbfgsbinp("lbfgsb.inp");
  //lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  //if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

  /*
  ConjugateGradientWSK CGsolver;//(true);
  int maxIter = 1000*model.dof();
  int restartStride = 2*model.dof();
  int printStride = 1000;//10*restartStride;
  double tol = 1.0e-6;//1.0e-8;
  double absTol = 1.0e-5; //1.0e-6;
  double tolLS = 1.0e-6;
  int maxIterLS = 20;
    
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 maxIter, restartStride, printStride, 
			 tol, absTol, tolLS, maxIterLS);
  CGsolver.setWolfeParameters(1.0e-6, 1.0e-3, 1.0e-4);

  //CGsolver.zeroOutData(true, false, false);
  */

  model.print("BeforeCG");
  
  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = bd.constraintVolume(); 
  double A = bd.constraintArea();
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
  std::cout << "Before CG:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Cons. Volume = "<< bd.constraintVolume() << std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Cons. Area = "<< bd.constraintArea() << std::endl
	    << "Reduced Volume = " << vred << std::endl
	    << "Cons. Reduced Volume = " << Vred << std::endl
	    << "Energy = " << solver.function() << std::endl
            <<  "ConstraintEnergy = " << bd.constraintEnergy() << std::endl
	    <<  "strainEnergy = " << bd.totalStrainEnergy() << std::endl;
    //<<  "workenergy = " << bd.workEnergy() <<std::endl;
   


 //model.checkRank(model.dof()-6,true);
  //  for(double nu=0.615424/*0.82127*/; nu > 0.6; nu-=0.1) {
//     bd.reduceVolume(nu); 
  for (double nu=0.98; nu>0.56; nu-=0.02){
  bd.reduceVolume(nu);
  //bd.resetReference();
  //model.checkConsistency(true, false);
  //return 0;    
  double energyRatio = 1.0;
  int iter=0, maxIters=10;
  //bd.printWork();
    
    while ( energyRatio > 1e-3 && iter < maxIters ) {	
      solver.solve( &model );      
      energyRatio = bd.constraintEnergy() / bd.totalStrainEnergy(); 
      iter++;
      bd.resetReference();
    }

        
    char name[20]; sprintf(name,"nu%f", nu);
    model.print(name);
    
    //bd.resetReference();
 
   
    v = bd.volume();
    a = bd.area(); 
    vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
    V = bd.constraintVolume();
    A = bd.constraintArea(); 
    Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
    
    std::cout << "nu = "<< nu << std::endl
      << "Volume = "<< bd.volume() << std::endl
      << "Cons. Volume = "<< bd.constraintVolume() << std::endl
      << "Area = "<< bd.area() << std::endl
      << "Cons. Area = "<< bd.constraintArea() << std::endl
      << "Reduced Volume = " << vred << std::endl
      << "Ref. Reduced Volume = " << Vred << std::endl
      << "Energy = " << solver.function() << std::endl
      << "ConstraintEnergy = " << bd.constraintEnergy() << std::endl
	      << "strainEnergy = " << bd.totalStrainEnergy() << std::endl;
      //<<  "workenergy = " << bd.workEnergy() <<std::endl;
    bd.printWork();
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
