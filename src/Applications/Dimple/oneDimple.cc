#include <string>
#include <iostream>
#include <vector>
#include <tvmet/Vector.h>
#include <fstream>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "TriangleQuadrature.h"
#include "LineElement.h"
#include "Spring.h"
#include "Model.h"
#include "Solver.h"
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

	  
  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;

  std::vector< int > dimpleNodesList;
  std::vector< int > dimpleElementsList;

  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;

    DeformationNode<3>::Point x;
    ifs >> id >> x(0) >> x(1) >> x(2);

    NodeBase::DofIndexMap idx(3);
    
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back(n);    
    defNodes.push_back(n); 
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

  int dimpleNodeNumber;
  int dimpleelementNumber;

  ifs >> dimpleNodeNumber >> dimpleelementNumber;
  for (int i=0; i<dimpleNodeNumber; i++){
    int nodeNumber;
    ifs >> nodeNumber;
    dimpleNodesList.push_back(nodeNumber-1);

  }
  
  for (int i=0; i<dimpleelementNumber; i++){
    int elementNumber;
    ifs >> elementNumber;
    dimpleElementsList.push_back(elementNumber-1);
  }

  ifs.close();


  //create the nodes and connectivities of dimple
  std::vector< DeformationNode<3>* > dimpleNodes;
  for(int i=0; i<dimpleNodesList.size(); i++){
    dimpleNodes.push_back(defNodes[dimpleNodesList[i]]);
  }

  std::vector< std::vector<int> > dimpleConnectivities;
  //vector< tvmet::Vector<int,3> > dimpleConnectivities;
  for(int i=0; i<dimpleElementsList.size(); i++) {
    std::vector<int> c(3);
    //tvmet::Vector<int, 3> c;
    for(int j=0; j<3; j++) c[j]=connectivities[dimpleElementsList[i]](j);
    dimpleConnectivities.push_back(c);
  }


  //create the node connnectivities for spring elements
  std::vector< tvmet::Vector<int,2> > springConnectivities;
  int scSize = 0;

  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> c(3);
    tvmet::Vector<int,2> n(2);

    for(int j=0; j<3; j++) c[j]=connectivities[i](j);
    for(int j=0; j<3; j++){

      bool pushFlag = true;
      n[0] = c[j];

      if (j < 2) 
	n[1] = c[j+1];
      else       
	n[1] = c[0];

      if (i==0) //don't need to check for the first element
	springConnectivities.push_back(n); 	
      else     //check if the spring already exists
	for(int k=0; k<scSize; k++){
	  if ((n[0] == springConnectivities[k][0] && n[1] == springConnectivities[k][1]) ||
	      (n[0] == springConnectivities[k][1] && n[1] == springConnectivities[k][0])   )
	    {
	      pushFlag = false; 
	      break;
	    }
	}
      if (pushFlag && i >0) springConnectivities.push_back(n);
    }
    scSize = springConnectivities.size();
  }


  double KC = 1.0e0;
  double KG = 0.0e0;
  double C0 = 0.0;
  double mu = 0.0;
  double KS = 0.0;
  double kSpring = 0.0;
  double Y = 1.0e2;
  double nu = 1.0/3.0;

  double dimpleSize = 0.1;        
  double epsilon = 0.0;
  double penaltyArea = 1.0e4;// penalty for dimple area constraint
  double chi = 20.0; // line tension
  double tau = 1.0; //surface tension 

  ifstream bdinp("constants.inp");
  bdinp >> dimpleSize >> epsilon >> kSpring >> penaltyArea >> Y >> chi >> tau;
  std::cout << "dimple size ="     << dimpleSize    << std::endl
	    << "epsilon ="     << epsilon    << std::endl
	    << "kSpring ="     << kSpring    << std::endl
	    << "penaltyArea =" << penaltyArea << std::endl
	    << "Y =" << Y << std::endl
	    << "Line tension=" << chi << std::endl
	    << "Surface tension=" << tau << std::endl;

  //pop up the dimple to locate on the sphere
  for(int i=0; i<dimpleNodeNumber; i++) {
    DeformationNode<3>::Point x;
    x = dimpleNodes[i]->point();

    double asquare = x(0)*x(0)+x(1)*x(1);

    double dimpleRadius = sqrt(dimpleSize*dimpleSize*(1.0+1.0/epsilon/epsilon));
    double b = sqrt(dimpleRadius*dimpleRadius - asquare);
    
    x(2)=b - dimpleSize/epsilon;
    
    dimpleNodes[i]->setPoint(x);
    dimpleNodes[i]->setPosition(x);
  }

  //find the nodes on the dimple boundary
  std::vector< int > lineNodesList;
  int lineNodeNumber;

  for(int i=0; i<dimpleNodeNumber; i++) {
    int n = dimpleNodesList[i];
    DeformationNode<3>::Point x;
    x = defNodes[n]->point();
    if (std::abs( norm2(x) - dimpleSize) < 1.0e-3)
      lineNodesList.push_back( n );
  }
  lineNodeNumber = lineNodesList.size();


  // create Body
  int quadOrder=2;

  // create Body for bending  
  EvansElastic bending( KC, KG, C0, mu, KS);
        
  typedef LoopShellBody<EvansElastic> LSB;
  //LSB bd(bending, connectivities, nodes, quadOrder,  
  // pressure, tension, totalCurvatureForce, penaltyVolume, penaltyArea, penaltyTotalCurvature);
  LSB bd(bending, connectivities, nodes, quadOrder,  
	 0.0, tau, 0.0, 0.0, 0.0, 0.0, noConstraint, multiplier, noConstraint);
  bd.setOutput(Body::paraview);

  // add springs to stabalize things
  std::vector<Spring<3>*> springs;
  for(int i=0; i<springConnectivities.size(); i++) {
    int m = springConnectivities[i](0);
    int n = springConnectivities[i](1);

    Spring<3> * sp = new Spring<3>(defNodes[m], defNodes[n], kSpring);
    springs.push_back( sp );
    bd.pushBack( sp );
  }

   std::cout << "Created a body for bending." << std::endl
	     << "  # nodes: " << bd.nodes().size()<<std::endl
	     << "  # elements: " << bd.elements().size() << std::endl
	     << "EvansElastic Material" <<std::endl
	     << " KG = " << KG << std::endl
	     << " C0 = " << C0 << std::endl
	     << " mu = " << mu << std::endl
	     << " KS = " << KS << std::endl
	     << " kSpring = " << kSpring << std::endl
	     << " surface tension = " << tau << std::endl;


  // create Body for stretching

  FVK stretching( 0.0, 0.0, C0, Y, nu );
  quadOrder = 2;
  typedef C0MembraneBody<TriangleQuadrature,FVK,ShapeTri3> MB; 
  MB bdm(stretching, dimpleConnectivities, nodes, quadOrder, 
	 0.0, 0.0, 0.0, penaltyArea, noConstraint, augmented);
  bdm.setOutput(Body::paraview);

  //create lineElement and pushback them to bdm
  for ( int i=0; i< lineNodeNumber; i++ ){
    int m = lineNodesList[i];
    int n = lineNodesList[i+1];
    if (i == lineNodeNumber - 1)
      n = lineNodesList[0];

    LineElement* le = new LineElement(chi, defNodes[m], defNodes[n]);

    bdm.pushBack(le);
  }

  //bdm.compute(true, true, false);  

  std::cout << "Created a body for stretching." << std::endl
	    << " # nodes: " << bdm.nodes().size()<<std::endl
	    << " # elements: " << bdm.elements().size() << std::endl
	    << " FVK Material " << std::endl
	    << " Y = " << Y << std::endl
	    << " nu = " << nu << std::endl
	    << " chi = " << chi << std::endl;
 

  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  
  bdc.push_back(&bdm);

  //
  // create Model
  Model model(bdc, nodes);

  //model.checkConsistency(true,false);
  //return 0;

  int m=5;
  double factr=0.0;//1.0e+7;
  double pgtol=1.0e-5;
  int iprint = -1;
  double pentol=1.0e-4;
  

  //ifstream lbfgsbinp("lbfgsb.inp");
  //lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  //if(verbose) 
  std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, 1000 );//(true);

//   double dt = 1.0e-8;
//   double tol = 1.0e-8;
//   double absTol = 1.0e-5;
//   int maxIteration = 10000;
//   int printStride = 2000;

//   ViscousRelaxation solver(dt, tol, absTol, maxIteration, printStride);

  // set bounds for solver
  bool boundaryConditionFlag = true;
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;
  if(boundaryConditionFlag) {
    //constrain the nodes at x, y directions ourside a certain distance
    for (int a=0; a<nodes.size(); a++){
      double x = nodes[a]->getPoint(0);
      double y = nodes[a]->getPoint(1);

      if (std::abs(x) > dimpleSize + 1.0 || std::abs(y) > dimpleSize + 1.0){
        for(int i=0; i<2; i++) {
          nbd(3*a+i) = 2;
          hi(3*a+i) = nodes[a]->getPoint(i);
          lo(3*a+i) = nodes[a]->getPoint(i);
        }
      }    
    }
  }

  
//   // set bounds for solver
//   bool boundaryConditionFlag = true;
//   blitz::Array<int,1> nbd(3*nodes.size());
//   blitz::Array<double,1> lo(3*nodes.size());
//   blitz::Array<double,1> hi(3*nodes.size());
//   nbd = 0;
//   lo = 0.0;
//   hi = 0.0;
//   if(boundaryConditionFlag) {
//     //constrain the nodes at x, y directions 
//     for (int a=0; a<nodes.size(); a++){
//       for(int i=0; i<2; i++) {
// 	nbd(3*a+i) = 2;
// 	hi(3*a+i) = nodes[a]->getPoint(i);
// 	lo(3*a+i) = nodes[a]->getPoint(i);
//       }
//     }
//     //constrain the dimple nodes
//     for(int a=0; a<dimpleNodeNumber; a++){
//        int b = dimpleNodesList[a];
//        nbd(3*b+2) = 2;
//        hi(3*b+2) = nodes[b]->getPoint(2);
//        lo(3*b+2) = nodes[b]->getPoint(2);
//     }
//   }

  
  solver.setBounds(nbd, lo, hi);

  double fT = 0.0;
//   bdm.updateFixedTension( fT );

  model.computeAndAssemble(solver,true,true,false);// initilize enregy and force
  
  double a = bd.area();
  double pda = bdm.prescribedArea();
  double da = bdm.area(); 
  std::cout << "Before CG:"<< std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Prescibed area=" <<bd.prescribedArea() <<std::endl
	    << "dimple area = " << da << std::endl
	    << "prescribed dimple are = " << pda << std::endl
	    << "bd0 Energy= " << bd.energy() << std::endl
    //<< "total length = " << bdm.totalLength() << std::endl
	    << "bd1 Energy = " << bdm.energy() << std::endl;

  //model.checkConsistency(true,false);
  //return 0;
  //bd.reduceVolume(0.98);


  double diffArea = 1.0;
  double energyRatio = 1.0;
  int iter=0, maxIter=10;
  
  model.print("BeforeCG");
  
  while (diffArea > 1.0e-3 || energyRatio > 1.0e-5  && iter <maxIter){
    bd.resetReference();
    solver.solve( &model );      

    double springEnergy = 0.0;
    for(int i=0; i<springs.size(); i++) springEnergy += springs[i]->energy();

    char name[30];sprintf(name,"oneDimple-%d", iter);
    model.print(name);
    
    a = bd.area();
    da = bdm.area(); 

    diffArea = std::abs(da - pda)/pda;
    energyRatio = springEnergy/bd.energy();

    std::cout  << "Area = "<< bd.area() << std::endl	       
	       << "dimple area = " << da << std::endl
	       << "area diff = " << diffArea << std::endl
	       << "bd0 Energy= " << bd.energy() << std::endl
	       << "energy ratio = " << energyRatio << std::endl
	       << "bd1 Energy = " << bdm.energy() << std::endl;


    iter++;

    fT += penaltyArea*(da-pda)/pda/pda;
    bdm.updateFixedTension( fT );

    if (diffArea > 1.0e-3){
      penaltyArea*=2.0;
      bdm.updatePenaltyArea( penaltyArea );
    }

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
