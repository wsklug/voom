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
#include "DimpleLoopShellBody.h"
#include "C0MembraneBody.h"
#include "C1LineElement.h"
#include "LineElement.h"
#include "ShapeTri3.h"
#include "TriangleQuadrature.h"
#include "ODQuadrature.h"
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

  std::vector< DeformationNode<3>* > nodesDimple;
  for(int i=0; i<dimpleNodesList.size(); i++){
    nodesDimple.push_back(defNodes[dimpleNodesList[i]]);
  }

  std::vector< std::vector<int> > connectivitiesDimple;
  //vector< tvmet::Vector<int,3> > connectivitiesDimple;
  for(int i=0; i<dimpleElementsList.size(); i++) {
    std::vector<int> c(3);
    //tvmet::Vector<int, 3> c;
    for(int j=0; j<3; j++) c[j]=connectivities[dimpleElementsList[i]](j);
    connectivitiesDimple.push_back(c);
  }

  std::vector< int > lineNodesList;
  int lineNodeNumber;

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
  bdinp >> dimpleSize >> epsilon >> kSpring>> penaltyArea >> Y >> chi >> tau;
  std::cout << "dimple size ="     << dimpleSize    << std::endl
	    << "epsilon ="     << epsilon    << std::endl
	    << "kSpring ="     << kSpring    << std::endl
	    << "penaltyArea =" << penaltyArea << std::endl
	    << "Y =" << Y << std::endl
	    << "Line tension =" << chi <<std::endl
	    << "surface tension=" << tau << std::endl;


  for(int i=0; i<dimpleNodeNumber; i++) {
    int n = dimpleNodesList[i];
    DeformationNode<3>::Point x;
    x = defNodes[n]->point();
    if (std::abs( x(0)*x(0)+x(1)*x(1) - dimpleSize*dimpleSize) < 1.0e-4)
      lineNodesList.push_back( n );
  }  
  lineNodeNumber = lineNodesList.size();


  //reorder the connectivities containing the line elelment
  std::vector< int > lineElementsList;
  for ( int i=0; i< lineNodeNumber; i++ ){
    int m = lineNodesList[i];
    int n = lineNodesList[i+1];
    if (i == lineNodeNumber - 1)
      n = lineNodesList[0];

    for(int a=0; a<connectivities.size(); a++) {
      std::vector<int> c(3);
      int third;
      bool change = false;

      for(int b=0; b<3; b++) c[b]=connectivities[a](b);



      for(int l=0; l<2; l++){ 
	for(int j=l+1; j<3; j++){

	  if ((c[l] == n && c[j] == m) || (c[l] == m && c[j] == n)){
	    for(int k=0; k<3; k++){
	      if (c[k] != m && c[k] != n) third = c[k];
	    }

	    DeformationNode<3>::Point x;
	    x = defNodes[third]->point();
	    if (tvmet::norm2(x) > dimpleSize) change = true;//choose the outside ring  
	  }
	}
      }

      if(change){

	connectivities[a][0] = n;
	connectivities[a][1] = m;
	connectivities[a][2] = third;

	lineElementsList.push_back(a);

	break;
      }

    }
  }



  //pop up the dimple
  for(int i=0; i<dimpleNodeNumber; i++) {
    DeformationNode<3>::Point x;
    x = nodesDimple[i]->point();

    double asquare = x(0)*x(0)+x(1)*x(1);

    double dimpleRadius = sqrt(dimpleSize*dimpleSize*(1.0+1.0/epsilon/epsilon));
    double b = sqrt(dimpleRadius*dimpleRadius - asquare);
    
    x(2)=b - dimpleSize/epsilon;
    
    nodesDimple[i]->setPoint(x);
    nodesDimple[i]->setPosition(x);
  }
  
  



  // create Body
  int nBoundaries=0;
  int quadOrder=2;

  // create Body for bending  
  EvansElastic bending( KC, KG, C0, mu, KS, kSpring);
        
  typedef LoopShellBody<EvansElastic> LSB;
  //LSB bd(bending, connectivities, nodes, quadOrder, nBoundaries, 
  // pressure, tension, totalCurvatureForce, penaltyVolume, penaltyArea, penaltyTotalCurvature, viscosity);
  LSB bd(bending, connectivities, nodes, quadOrder, nBoundaries, 
	 0.0, tau, 0.0, 0.0, 0.0, 0.0, 0.0);
  bd.setOutput(Body::paraview);


  //create C1LineElement and pushback them to bd
//   LSB::FeElementContainer shells = bd.shells();
//   ODQuadrature quad1(2);
//   for(int i=0; i<lineElementsList.size(); i++){
//     C1LineElement<EvansElastic>* le = new C1LineElement<EvansElastic>(shells[lineElementsList[i]], quad1, chi);
//     bd.pushBackC1LineElement(le);
//   }


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
  MB bdm(stretching, connectivitiesDimple, nodes, quadOrder, 0.0, 0.0, penaltyArea);
  bdm.setOutput(Body::paraview);

  //create lineElement and pushback them to bdm
  for ( int i=0; i< lineNodeNumber; i++ ){
    int m = lineNodesList[i];
    int n = lineNodesList[i+1];
    if (i == lineNodeNumber - 1)
      n = lineNodesList[0];

    LineElement* le = new LineElement(chi, defNodes[m], defNodes[n]);

    bdm.pushBackLineElement(le);

  }

  //bdm.compute(true, true, false);  

  //std::cout << "length =" << bdm.totalLength()<< std::endl;
  //bdm.setMaterialY (10.0);
  std::cout << "Created a body for stretching." << std::endl
	    << " # nodes: " << bdm.nodes().size()<<std::endl
	    << " # elements: " << bdm.elements().size() << std::endl
	    << " FVK Material " << std::endl
	    << " Y = " << Y << std::endl
	    << " nu = " << nu << std::endl
	    << " chi = " << chi << std::endl;
  
  //bdm.setConstraintArea(1.0);

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
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

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

      if (std::abs(x) > dimpleSize + 2.0 || std::abs(y) > dimpleSize + 2.0){ 
	for(int i=0; i<2; i++) {
	  nbd(3*a+i) = 2;
	  hi(3*a+i) = nodes[a]->getPoint(i);
	  lo(3*a+i) = nodes[a]->getPoint(i);
	}
      }    
    }
  }

//   if(boundaryConditionFlag) {
//     //constrain all the nodes to move only in z direction
//     for (int a=0; a<nodes.size(); a++){
//       for(int i=0; i<2; i++) {
// 	nbd(3*a+i) = 2;
// 	hi(3*a+i) = nodes[a]->getPoint(i);
// 	lo(3*a+i) = nodes[a]->getPoint(i);
//       }    
//     }

//     //no constraint on the dimple nodes
//     for(int a=0; a<dimpleNodeNumber; a++){
//        int b = dimpleNodesList[a];
//        nbd(3*b) = 0;
//        nbd(3*b+1) = 0;
//        nbd(3*b+2) = 0;
//     }
//   }

  
  solver.setBounds(nbd, lo, hi);

  double fT = 0.0;
  bdm.updateFixedTension( fT );

  model.computeAndAssemble(solver,true,true,false);// initilize enregy and force
  
  double a = bd.area();
  double cda = bdm.constraintArea();
  double da = bdm.area(); 
  std::cout << "Before CG:"<< std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Constraint area=" <<bd.constraintArea() <<std::endl
	    << "dimple area = " << da << std::endl
	    << "constraint dimple are = " << cda << std::endl
	    << "bd0 Energy= " << bd.energy() << std::endl
	    << "total length = " << bdm.totalLength() << std::endl
	    << "bd1 Energy = " << bdm.energy() << std::endl;

  //model.checkConsistency(true,false);
  //return 0;
  //bd.reduceVolume(0.98);


  double diffArea = 1.0;
  double energyRatio = 1.0;
  //double energyRatio2 = 1.0;
  int iter=0, maxIter=10;
  //double en=0;
  //double enplus; 
  
  //  model.print("tension");  
  
  
  model.print("BeforeCG1");
  
  while (/*diffArea > 1.0e-3 || energyRatio > 1.0e-4 && bd.constraintEnergy() > 1.0e-5 &&*/ iter <maxIter){
    bd.resetReference();


    solver.solve( &model );      

    char name[30];sprintf(name,"oneDimple-%d", iter);
    model.print(name);
    
    a = bd.area();
    da = bdm.area(); 

    diffArea = std::abs(da - cda)/cda;
    energyRatio = bd.viscousEnergy()/bd.energy();
    energyRatio = std::abs(energyRatio);
    //energyRatio2 = bdm.totalStrainEnergy()/bdm.energy();

    std::cout  << "Area = "<< bd.area() << std::endl	       
	       << "dimple area = " << da << std::endl
	       << "area diff = " << diffArea << std::endl
	       << "bd0 Energy= " << bd.energy() << std::endl
	       << "bd0 constraint Energy=" << bd.constraintEnergy()<<std::endl
	       << "bd0 line Energy=" << bd.lineEnergy()<<std::endl
	       << "energy ratio = " << energyRatio << std::endl
	       << "bd1 Energy = " << bdm.energy() << std::endl;
// 	       << "bd1 constraint Energy = " << bdm.totalStrainEnergy() << std::endl
// 	       << "energy ratio 2= " << energyRatio2 << std::endl;
	       
    //enplus = bd.energy();
    //energyRatio = (enplus-en)/enplus;
    //en = enplus;

    iter++;


    fT += penaltyArea*(da-cda)/cda/cda;
    bdm.updateFixedTension( fT );

    std::cout << "tension force of dimple = " << fT << std::endl;

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
