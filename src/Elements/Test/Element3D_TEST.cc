// Test for birck elements written by Luigi Perotti from
// Melissa M. Gibbons test for NonlinearElastic element class

#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include "Element3D.h"
#include "NonlinearElastic.h"
#include "CompNeoHookean.h"
#include "StVenant.h"
#include <iostream>

using namespace std;
using namespace voom;

int main()
{ 
  // Make a 4 nodes tet element
  std::vector< DeformationNode<3>* > Nodes;
  Nodes.reserve(4);
  std::vector< DeformationNode<3>* > nNodes;
  nNodes.reserve(4);


  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 0.0; 
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( a );
  nNodes.push_back( a );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 1.0, 0.0; 
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( b );
  nNodes.push_back( b ); 

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 0.0; 
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( c );
  nNodes.push_back( c );


  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 1.0; 
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( d );
  nNodes.push_back( d );

  // create material
  double rho = 500.0;
  double E = 200.0;
  double nu = 0.4;
  CompNeoHookean protein( rho, E, nu );



  // create element
  TetQuadrature Quad(1);
  ShapeTet4 Sh;
 
  d->setPoint(2,1.1);
  Element3D Tet(Nodes, &protein, Quad, Sh);
  /*cout << "Node a force " <<  a->force() << endl;
  cout << "Node b force " <<  b->force() << endl;
  cout << "Node c force " <<  c->force() << endl;
  cout << "Node d force " <<  d->force() << endl;
  cout << "energy = " << Tet.energy() << endl;
  */
  cout << "I initialized the element ok" << endl;

  NonlinearElastic<TetQuadrature, CompNeoHookean, ShapeTet4> RefEl(Quad, protein, nNodes);

  Tet.checkConsistency();
  RefEl.checkConsistency();

  cout << endl << "Invariants: " << endl;
  int f = 0;
  vector<pair<Vector3D, vector<double > > > Invariants = Tet.invariants(f);
  for (int i = 0; i<Invariants.size(); i++) {
    cout << (Invariants[i]).first << "  ";
    for (int j =0; j<3; j++) {
      cout << ((Invariants[i]).second)[j] << " "; 
    }
    cout << endl;
  }
  cout << endl;

  return 0;
}
