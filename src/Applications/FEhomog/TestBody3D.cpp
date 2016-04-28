#include <string>
#include <iostream>
#include <vector>
#include "Body3D.h"
#include "Capsid3DBody.h"
#include "Node.h"
#include "CompNeoHookean.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
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
  CompNeoHookean protein1( rho, E, nu );
  CompNeoHookean protein2( rho, E, nu );

  vector<Material *> VecProtein;
  VecProtein.push_back(&protein1);
  VecProtein.push_back(&protein2);



  // Connectivity table
  vector<vector< int> > Connectivities(2, vector<int>(4, 0));
  Connectivities[0][0] = 0; Connectivities[0][1] = 1; Connectivities[0][2] = 2; Connectivities[0][3] = 3;
  Connectivities[1][0] = 0; Connectivities[1][1] = 1; Connectivities[1][2] = 2; Connectivities[1][3] = 4;



  // DefNodeContainer
  vector<DeformationNode<3> *> Nodes;
  Nodes.reserve(5);
  Body::NodeContainer CapsidNodes;
  CapsidNodes.reserve(5);

  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 0.0; 
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( a );
  CapsidNodes.push_back( a );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 1.0, 0.0; 
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( b );
  CapsidNodes.push_back( b );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 0.0; 
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( c );
  CapsidNodes.push_back( c );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 1.0; 
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( d );
  CapsidNodes.push_back( d );
  
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0,-1.0; 
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( e );
  CapsidNodes.push_back( e );
 

  // Create quadrature rule and shape function
  TetQuadrature Quad(1);
  ShapeTet4 Sh;
 
  d->setPoint(2,1.1);
  e->setPoint(1,-0.6);
  
  // Create body
  Body3D BodyTest(VecProtein, Connectivities, Nodes, Quad, Sh);
  /*cout << "I initialized the body" << endl;
  cout << "Node a force " <<  a->force() << endl;
  cout << "Node b force " <<  b->force() << endl;
  cout << "Node c force " <<  c->force() << endl;
  cout << "Node d force " <<  d->force() << endl;*/

 
  typedef Capsid3DBody<TetQuadrature, CompNeoHookean, ShapeTet4> Capsid;
  Capsid CapsidBody(protein1, Connectivities, CapsidNodes, 2);
  /*cout << "Node a force " <<  a->force() << endl;
  cout << "Node b force " <<  b->force() << endl;
  cout << "Node c force " <<  c->force() << endl;
  cout << "Node d force " <<  d->force() << endl;*/

  /*for(Body::NodeIterator n = CapsidNodes.begin(); n != CapsidNodes.end(); n++) 
	for(int i=0; i<(*n)->dof(); i++)
	(*n)->setForce(i,0.0);*/
  
  BodyTest.checkConsistency();
  CapsidBody.checkConsistency();
  


  cout << endl << "Invariants: " << endl;
  int f = 0;
  vector<pair<Vector3D, vector<double > > > Invariants = BodyTest.invariants(f);
  for (int i = 0; i<Invariants.size(); i++) {
    cout <<  Invariants[i].first << endl;
    for (int j =0; j<3; j++) {
      cout <<  (Invariants[i].second)[j] << " "; 
    }
    cout << endl;
  }
  cout << endl;

  for (int k = 0; k<Nodes.size(); k++)
  {
    delete Nodes[k];
  }
  
	
  return 0;
}



