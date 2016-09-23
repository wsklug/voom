#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "StVenant.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "Capsid3DBody.h"


// #define NODENUMBER 42
// #define ELEMNUMBER 80

using namespace tvmet;
using namespace std;
using namespace voom;

//void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{

  ifstream ifs;
  string ofn;
 
  string modelName = argv[1];
  string inputFileName = modelName + ".vtk";

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
  char key;  ifs>>key;
  
  // MMG, 9-18-06:  Input .vtk file containing nodes and connectivities
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  ifs >> token;// skip number type

  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
  }
  cout << "Number of nodes: " <<nodes.size() << endl;

  // read in tetrahedral connectivities
  while( token != "CELLS" ) ifs >> token;
  vector< tvmet::Vector<int,4> > connectivities;
  tvmet::Vector<int, 4> c;
  int ntet=0; ifs >> ntet;
  connectivities.reserve(ntet);
  cout << "Number of tetrahedrons: " << ntet << endl;
  int ntmp=0;
  ifs >> ntmp;
  if(ntmp != 5*ntet) {
    cout << "Non-tetrahedral elements?" << endl;
    return 0;
  }
  for (int i = 0; i<ntet; i++){
    int tmp=0;
    ifs >> tmp;
    ifs >> tmp; c[0]=tmp;
    ifs >> tmp; c[1]=tmp;
    ifs >> tmp; c[2]=tmp;
    ifs >> tmp; c[3]=tmp;
    connectivities.push_back(c);
  }
  if(verbose) cout << connectivities.size() << endl;
  ifs.close();

  // create material
  typedef StVenant MaterialType;
  MaterialType protein( 500.0, 200.0, 0.4 );
  std::cout << "Material has been created." << std::endl;
 
  //
  // create Body
  typedef Capsid3DBody<MaterialType> Capsid;
  Capsid bd(protein, connectivities, nodes, 2);

  bd.setOutput(Body::paraview);
  cout << "Created a body." << endl;
 
  bd.compute(false,true, false);
  bd.printParaview("staticForce");
 
  
  // bd.rankTest();
  //	bd.checkConsistency();
	
  //	ifs.close();
	
    return 0;
}


// void ioSetting(int argc, char* argv[], ifstream& ifs, string& ofn)
// {

//   // if no arguement or too many
//   if (argc < 2 || argc > 3)
//     {
//       cout << "Usage: ProgramName InputFileName [OutputFileName]." << endl;
//       exit(0);
//     }

//   string inFullName = argv[1];
//   string pathName;
//   //	basic_string <char>::size_type sztp = inFullName.find_last_of("/");
	
//   // 	if ( (sztp) != string::npos )  // no position was found
//   // 		pathName = inFullName.substr(0, sztp) + "/";
//   // 	else
//   // 		pathName = "";

//   pathName = "";
	
//   // only one arguement
//   if ( argc == 2 ) ofn = pathName + "output.dat";


//   // two arguements
//   if ( argc == 3 ) {
//     ofn = argv[2];
//     ofn = pathName + ofn;
//   }

//   // create input stream
//   ifs.open( inFullName.c_str(), ios::in);
//   if (!ifs)
//     {
//       cout << "can not open input file: " << inFullName << endl;
//       exit(0);
//     }
	
// }
