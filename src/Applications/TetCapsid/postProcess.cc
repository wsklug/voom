// -*- C++ -*-
//----------------------------------------------------------------------
//
//                Melissa M. Gibbons & William S. Klug
//                University of California Los Angeles
//                   (C) 2007 All Rights Reserved
//
//----------------------------------------------------------------------

// Program that will take in an indentation output file with: nodes,
// elements, and displacements, and will calculate stresses, strains,
// and internal forces.

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
#include "ShapeTet4CP.h"
#include "Model.h"

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{

  bool verbose=true;

  string modelName = argv[1]; 
  string inputFileName = modelName + ".vtk";

  // create input stream
  ifstream ifs;
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
  
  // Input .vtk file containing nodes and connectivities
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
  std::vector<int> c(4);
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
  
  while( token != "displacements" ) ifs >> token;
  ifs >> token;// skip number type
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    x += defNodes[i]->position();
    defNodes[i]->setPoint(x);
  }
  
  ifs.close();
  
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
  unsigned int quadOrder = 1;
  
  typedef Capsid3DBody<TetQuadrature,MaterialType,ShapeTet4> Capsid;
  Capsid bd(protein, connectivities, nodes, quadOrder);

  bd.setOutput(Body::paraview);

  string name = modelName + "PostProcess";
  bd.PostProcess( name );

  std::cout << "All done.  Bye now." << std::endl;

  return 0;
}

