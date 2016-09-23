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
#include "SpringBody.h"
#include "ShapeTri3.h"
#include "TriangleQuadrature.h"
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
//   for(int i=0; i<dimpleNodesList.size(); i++)
//   cout<< dimpleNodesList[i]<< endl;
  for (int i=0; i<dimpleelementNumber; i++){
    int elementNumber;
    ifs >> elementNumber;
    dimpleElementsList.push_back(elementNumber-1);
  }

  ifs.close();



  double dimpleSize = 0.1;
  double epsilon = 0.15;
  double center1[2]={-0.2,0};
  double center2[2]={0.2,0};
  double tilt1 = 0.0;
  double tilt2 = 0.0;

  double KC = 1.0e0;
  double KG = 0.0e0;
  double C0 = 0.0;
  double mu = 0.0;
  double KS = 0.0;
  double kSpring = 0.0;
  double viscosity = 0.0e0;
  //double Y = 5.0e5;
  //double nu = 1.0/3.0;
  double springConstant = 1.0e3;

  ifstream bdinp("constants2D.inp");
  bdinp >> dimpleSize >> epsilon >> center1[0] >> center2[0] >> springConstant;
  std::cout << "dimple size ="     << dimpleSize    << std::endl
	    << "epsilon = "        << epsilon       << std::endl
	    << "dimple 1 center: " << center1[0] << "," << 0 << std::endl        
	    << "dimple 2 center: " << center2[0] << "," << 0 << std::endl
	    << "SpringConstant k = " << springConstant << std::endl;  


  //create the node list for spring element
  std::vector< tvmet::Vector<int,2> > SpringConnectivities;
  int scSize = 0;

  for(int i=0; i<dimpleElementsList.size(); i++) {
    std::vector<int> c(3);
    tvmet::Vector<int,2> n(2);

    for(int j=0; j<3; j++) c[j]=connectivities[dimpleElementsList[i]](j);
    for(int j=0; j<3; j++){

      bool pushFlag = true;
      n[0] = c[j];

      if (j < 2) 
	n[1] = c[j+1];
      else       
	n[1] = c[0];


      if (i==0) 
	SpringConnectivities.push_back(n); 	

      else
	for(int k=0; k<scSize; k++){
	  if ((n[0] == SpringConnectivities[k][0] && n[1] == SpringConnectivities[k][1]) ||
	      (n[0] == SpringConnectivities[k][1] && n[1] == SpringConnectivities[k][0])   )
	    {
	      pushFlag = false; 
	      break;
	    }
	}
      if (pushFlag && i >0) SpringConnectivities.push_back(n);

    }

    scSize = SpringConnectivities.size();
  }


  std::vector< DeformationNode<3>* > nodesDimple;
  for(int i=0; i<dimpleNodesList.size(); i++){
    nodesDimple.push_back(defNodes[dimpleNodesList[i]]);
  }
  
  //pop up the two dimple
  for(int i=0; i<dimpleNodeNumber; i++) {
    DeformationNode<3>::Point x;
    x = nodesDimple[i]->point();

    double asquare1 = (x(0)-center1[0])*(x(0)-center1[0])+(x(1)-center1[1])*(x(1)-center1[1]);
    double asquare2 = (x(0)-center2[0])*(x(0)-center2[0])+(x(1)-center2[1])*(x(1)-center2[1]);

    double asquare = asquare1;
    if (asquare1 > asquare2) asquare = asquare2;

    double dimpleRadius = sqrt(dimpleSize*dimpleSize*(1.0+1.0/epsilon/epsilon));
    double b = sqrt(dimpleRadius*dimpleRadius - asquare);
    
    x(2)=b - dimpleSize/epsilon;
    
    nodesDimple[i]->setPoint(x);
    nodesDimple[i]->setPosition(x);
  }

    

  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=1.0;
  double totalCurvatureForce=0.0;
  int quadOrder=1;
  double penaltyVolume=0.0e4;
  double penaltyArea=0.0e4;
  double penaltyTotalCurvature=0.0e4;

  // create Body for bending  
  EvansElastic bending( KC, KG, C0, mu, KS, kSpring);
        
  typedef LoopShellBody<EvansElastic> LSB;
  LSB bd(bending, connectivities, nodes, quadOrder, nBoundaries, 
	 pressure, tension, totalCurvatureForce, penaltyVolume, penaltyArea, penaltyTotalCurvature, viscosity); 
   bd.setOutput(Body::paraview);



   std::cout << "Created a body for bending." << std::endl
	     << "  # nodes: " << bd.nodes().size()<<std::endl
	     << "  # elements: " << bd.elements().size() << std::endl
	     << "EvansElastic Material" <<std::endl
	     << " KG = " << KG << std::endl
	     << " C0 = " << C0 << std::endl
	     << " mu = " << mu << std::endl
	     << " KS = " << KS << std::endl
	     << " kSpring = " << kSpring << std::endl
	     << " tension = " << tension << std::endl;

  

  // create Body for stretching
  typedef SpringBody MB; 
  MB bdm(nodes, SpringConnectivities, springConstant);
  bdm.setOutput(Body::paraview);

  std::cout << "Created a body for springs." << std::endl;

  

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

  //get four nodes at two dimples to compute the tilting angle
  int fourNodes[4];
  for(int a=0; a<dimpleNodeNumber; a++){
    int b = dimpleNodesList[a];
    double x = nodes[b]->getPoint(0);
    double y = nodes[b]->getPoint(1);
    double distance1 = sqrt((x-(center1[0]+dimpleSize))*(x-(center1[0]+dimpleSize)) + (y-center1[1])*(y-center1[1]));
    double distance2 = sqrt((x-(center1[0]-dimpleSize))*(x-(center1[0]-dimpleSize)) + (y-center1[1])*(y-center1[1]));
    double distance3 = sqrt((x-(center2[0]-dimpleSize))*(x-(center2[0]-dimpleSize)) + (y-center2[1])*(y-center2[1]));
    double distance4 = sqrt((x-(center2[0]+dimpleSize))*(x-(center2[0]+dimpleSize)) + (y-center2[1])*(y-center2[1]));

    if(distance1 < 0.001) fourNodes[0]=b;
    if(distance2 < 0.001) fourNodes[1]=b;
    if(distance3 < 0.001) fourNodes[2]=b;
    if(distance4 < 0.001) fourNodes[3]=b;
  }

  std::cout << "Four nodes on two dimples are " << std::endl
	    << "Node " << fourNodes[0] << " at x = " << nodes[fourNodes[0]]->getPoint(0) << " y = "  << nodes[fourNodes[0]]->getPoint(1) << std::endl
	    << "Node " << fourNodes[1] << " at x = " << nodes[fourNodes[1]]->getPoint(0) << " y = "  << nodes[fourNodes[1]]->getPoint(1) << std::endl
            << "Node " << fourNodes[2] << " at x = " << nodes[fourNodes[2]]->getPoint(0) << " y = "  << nodes[fourNodes[2]]->getPoint(1) << std::endl
	    << "Node " << fourNodes[3] << " at x = " << nodes[fourNodes[3]]->getPoint(0) << " y = "  << nodes[fourNodes[3]]->getPoint(1) << std::endl;

  // set bounds for solver
  bool boundaryConditionFlag = true;
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;
  if(boundaryConditionFlag) {
    //constrain all the nodes to move only in z direction
    for (int a=0; a<nodes.size(); a++){
      for(int i=0; i<2; i++) {
	nbd(3*a+i) = 2;
	hi(3*a+i) = nodes[a]->getPoint(i);
	lo(3*a+i) = nodes[a]->getPoint(i);
      }    
    }
   
    //now let nodes on dimples to move freely between center2 + dimpleSize < x < center2 + dimpleSize 
    //except two closet one
    for(int a=0; a<dimpleNodeNumber; a++){
       int b = dimpleNodesList[a];
       nbd(3*b) = 2;
       nbd(3*b+1) = 2;
       hi(3*b) = center2[0] + dimpleSize;
       lo(3*b) = center1[0] - dimpleSize;
       hi(3*b+1) = nodes[b]->getPoint(1);
       lo(3*b+1) = nodes[b]->getPoint(1);

       nbd(3*b+2) = 0;
       

       if (b == fourNodes[0] || b == fourNodes[2]){
	 for(int i=0; i<3; i++) {	  
	   nbd(3*b+i) = 2;
	   hi(3*b+i) = nodes[b]->getPoint(i);
	   lo(3*b+i) = nodes[b]->getPoint(i);
	 }   
	 cout << "fixing node " << b << " at x = " << nodes[b]->getPoint(0) << " y = " << nodes[b]->getPoint(1) << endl;
       }


    }
  }
  
  solver.setBounds(nbd, lo, hi);

  model.print("BeforeCG");

  //exit(0);
  double a = bd.area(); 
  std::cout << "Before CG:"<< std::endl
	    << "Area = "<< bd.area() << std::endl;

  //model.checkConsistency(true,false);
  //return 0;
  //bd.reduceVolume(0.98);
  double energyRatio = 1.0;
  int iter=0, maxIter=10;
  double en=0;
  double enplus; 

  while (energyRatio > 10e-6 && iter <maxIter){
    solver.solve( &model );      

    char name[30];sprintf(name,"twoDimple");
    model.print(name);
    
    a = bd.area(); 

    std::cout  << "Area = "<< bd.area() << std::endl
	       << "bd0 Energy= " << bd.energy() << std::endl
	       << "bd1 Energy = " << bdm.energy() << std::endl;


    enplus = bd.energy();
    energyRatio = (enplus-en)/enplus;
    en = enplus;

    vector< DeformationNode<3>::Point > fourPoints;
    for (int i=0; i<4; i++){
      int a = fourNodes[i];
      DeformationNode<3>::Point  points;
      points[0] = nodes[a]->getPoint(0);
      points[1] = nodes[a]->getPoint(1);
      points[2] = nodes[a]->getPoint(2);
      fourPoints.push_back(points);
    }
    double dz = fourPoints[0][2] - fourPoints[1][2];
    double xy = sqrt((fourPoints[0][0] - fourPoints[1][0])*(fourPoints[0][0] - fourPoints[1][0]) +
		     (fourPoints[0][1] - fourPoints[1][1])*(fourPoints[0][1] - fourPoints[1][1]));
    tilt1 = dz/xy;

    dz = fourPoints[2][2] - fourPoints[3][2];
    xy = sqrt((fourPoints[2][0] - fourPoints[3][0])*(fourPoints[2][0] - fourPoints[3][0]) +
	      (fourPoints[2][1] - fourPoints[3][1])*(fourPoints[2][1] - fourPoints[3][1]));
    tilt2 = dz/xy;

    std::cout << "tilting anlge alpha = " << tilt1 << std::endl
	      << "two nodes are: Node " << fourNodes[0] << " at " << fourPoints[0][0] <<"," << fourPoints[0][1] << "," << fourPoints[0][2]
              << " Node " << fourNodes[1] << " at " << fourPoints[1][0] <<"," << fourPoints[1][1] << "," << fourPoints[1][2] << std::endl
	      << "tilting anlge beta = " << tilt2 << std::endl
	      << "two nodes are: Node " << fourNodes[2] << " at " << fourPoints[2][0] <<"," << fourPoints[2][1] << "," << fourPoints[2][2]
              << " Node " << fourNodes[3] << " at " << fourPoints[3][0] <<"," << fourPoints[3][1] << "," << fourPoints[3][2] << std::endl;


//     for(int a=0; a<dimpleNodeNumber; a++){
//        int b = dimpleNodesList[a];
//        nbd(3*b) = 2;
//        nbd(3*b+1) = 2;
//        hi(3*b) =  nodes[fourNodes[2]]->getPoint(0) + 2.0*dimpleSize*cos(tilt2/1.4);
//        lo(3*b) =  nodes[fourNodes[0]]->getPoint(0) - 2.0*dimpleSize*cos(tilt1/1.4);
//        hi(3*b+1) = nodes[b]->getPoint(1);
//        lo(3*b+1) = nodes[b]->getPoint(1);

//        nbd(3*b+2) = 0;
       

//        if (b == fourNodes[0] || b == fourNodes[2]){
// 	 for(int i=0; i<3; i++) {	  
// 	   nbd(3*b+i) = 2;
// 	   hi(3*b+i) = nodes[b]->getPoint(i);
// 	   lo(3*b+i) = nodes[b]->getPoint(i);
// 	 }   
//        }

//     }

//     solver.setBounds(nbd, lo, hi);


    iter++;

    bdm.updateFixedForce();
    bdm.updateSpringConstatnt(3.0, 1.0e-4);

    string egName =  "twoDimple.eg";
    ofstream eg(egName.c_str());

    double distance = nodes[fourNodes[2]]->getPoint(0) - nodes[fourNodes[0]]->getPoint(0) + 2.0*dimpleSize; 
    eg   << std::setw( 16 ) << std::setprecision(12) <<  distance
	 << std::setw( 16 ) << std::setprecision(12) <<  bd.energy()
	 << std::setw( 16 ) << std::setprecision(12) <<  tilt1
	 << std::setw( 16 ) << std::setprecision(12) <<  tilt2
	 << std::endl;

    eg.close();

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
