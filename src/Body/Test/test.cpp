#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "SCElastic.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"


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

	
  int ngp = 3;
  // 	cout << "input number of gauss point(1 or 3 ONLY): " << endl;
  // 	cin >> ngp;
  if ( ngp != 1 && ngp != 3){
    cout << "incorrect # of gauss point: input 1 or 3 please." << endl;
    exit(0);
  }

  int NODENUMBER, ELEMNUMBER;
	
  ifs >> NODENUMBER >> ELEMNUMBER;
	
  //
  // create vector of nodes
  std::vector< NodeBase* > nodes;
  for ( int i = 0; i < NODENUMBER; i++){
    int id=-1;
    DeformationNode<3>::Point x;
    ifs >> id >> x(0) >> x(1) >> x(2);
    cout << setw(12) << id 
	 << setw(20) << x(0)
	 << setw(20) << x(1)
	 << setw(20) << x(2) << endl;
    NodeBase::DofIndexMap idx(3);
    nodes.push_back(new DeformationNode<3>(id,idx,x));
  }
  cout << nodes.size() << endl;
  for(int i=0; i<NODENUMBER; i++) {
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
	
  //                           KG
  //      KC___________         |       ________________Spontaneous curvature
  //                   |        |       |
  //                   |        |       |
  //                   |        |       |
  //                   |        |       |
  //                   |        |       |
  //                   |        |       |
  //                   V        V       V
  SCElastic bilayer( 1.0e02, 0.0e00, 5.0e00 );
  std::cout << "SCElastic Material has been created." << std::endl;
	
  //
  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double penaltyVolume=0.0;
  double penaltyArea=0.0;
  LoopShellBody<SCElastic> bd(bilayer, connectivities, nodes, ngp, 
			      nBoundaries, pressure, penaltyVolume, penaltyArea);
  //
  // 
  cout << "Created a body." << endl;
  // bd.printByHDS();
  bd.compute(false,true, false);
  bd.printParaview("staticForce");
  //	bd.createOpenDXData(ofn, 0);
  //	bd.createInputFile("InputFile");
  cout << "volume of current body = " << bd.volume() << endl;
  //	cout << " rank test ..." << endl;
  
  // bd.rankTest();
  //	bd.checkConsistency();
  //	bd.ComputeVescileForce();
	
    //	ifs.close();
	
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
