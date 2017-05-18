#include <string>
#include <iostream>
#include <vector>
#include "ViscosityBody.h"

#include "Node.h"

#include <tvmet/Vector.h>
#include <fstream>



using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  // Create nodes
  vector<DeformationNode<3> *> Nodes;
  Nodes.reserve(5);

  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point X, x;
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  X = 0.0, 0.0, 0.0; 
  x = X*1.1;
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, X, x);
  Nodes.push_back( a );
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  X = 0.0, 1.0, 0.0; 
  x = X*1.1;
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, X, x);
  Nodes.push_back( b );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  X = 1.0, 0.0, 0.0; 
  x = X*1.1;
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, X, x);
  Nodes.push_back( c );
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  X = 0.0, 0.0, 1.0; 
  x = X*1.1;
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, X, x);
  Nodes.push_back( d );
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  X = 0.0, 0.0,-1.0; 
  x = X*1.1;
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, X, x);
  Nodes.push_back( e );

  vector<tvmet::Vector<int, 3> > Connectivity(3, tvmet::Vector<int, 3>(0));
  Connectivity[0](0) = 0; Connectivity[0](1) = 1; Connectivity[0](2) = 2;
  Connectivity[1](0) = 0; Connectivity[1](1) = 1; Connectivity[1](2) = 3;
  Connectivity[2](0) = 0; Connectivity[2](1) = 3; Connectivity[2](2) = 4;
  Connectivity[2](0) = 1; Connectivity[2](1) = 3; Connectivity[2](2) = 4;

  // Spring constants
  const double sigma   = 1.76;

  // Create body
  ViscosityBody TestBody(Nodes, Connectivity, sigma);

  //Perturb nodal positions
  x = 0.5, 0.5, 0.5;
  a->setPoint(x);
  x = 1.5, 0.5, -0.52;
  c->setPoint(x);

  // Check consistency
  TestBody.checkConsistency(false); // verbose set to false, for more details set it to true

  // Print Body Energy
  TestBody.compute(true, false, false);
  cout << "Viscosity body energy = " << TestBody.totalStrainEnergy() << endl;

  TestBody.resetRefConf();
  TestBody.compute(true, false, false);
  cout << "Viscosity body energy after reset = " << TestBody.totalStrainEnergy() << endl;
 	
  return 0;
}



