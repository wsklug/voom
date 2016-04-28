#include <string>
#include <iostream>
#include <vector>
#include "Body3D.h"
//#include "Capsid3DBody.h"
#include "Node.h"
#include "CompNeoHookean.h"
#include "TetQuadrature.h"
#include "ShapeTet4.h"
#include <tvmet/Vector.h>
#include <fstream>


using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{

  // Create material
  double rho = 500.0;
  double E = 200.0;
  double nu = 0.4;
  CompNeoHookean protein( rho, E, nu );



  // Connectivity table
  vector<vector< int> > Connectivities(2, vector<int>(4, 0));
  Connectivities[0][0] = 0; Connectivities[0][1] = 1; Connectivities[0][2] = 2; Connectivities[0][3] = 3;
  Connectivities[1][0] = 0; Connectivities[1][1] = 1; Connectivities[1][2] = 2; Connectivities[1][3] = 5;



  // DefNodeContainer
  vector<DeformationNode<3> *> Nodes;
  Nodes.reserve(5);

  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 0.0; 
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( a );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 1.0, 0.0; 
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( b );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 0.0; 
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( c );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 1.0; 
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( d );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0,-1.0; 
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( e );

 

  // Create quadrature rule and shape function
  TetQuadrature Quad(1);
  ShapeTet4 Sh;
 
  
  
  // Create body
  Body3D BodyTest(&protein, Connectivities, Nodes, Quad, Sh);
  cout << "I initialized the body" << endl;
 
  // typedef Capsid3DBody<TetQuadrature, CompNeoHookean, ShapeTet4> Capsid;
  // Capsid bd(protein, Connectivities, Nodes, 2);
  
  d->setPoint(2,1.1);
  BodyTest.checkConsistency();
  /*
  cout << endl << "Invariants: " << endl;
  vector<vector<double > > Invariants = Tet.invariants();
  for (int i = 0; i<Invariants.size(); i++) {
    for (int j =0; j<3; j++) {
      cout << Invariants[i][j] << " "; 
    }
    cout << endl;
  }
  cout << endl;
  */
	
  return 0;
}



