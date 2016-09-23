#include <string>
#include <iostream>
#include <vector>


#include "Node.h"
#include "PotentialBody.h"
#include "ProteinLennardJones.h"
#include "ProteinBody.h"

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

  vector<ProteinNode *> ProteinNodes;
  ProteinNodes.reserve(5);

  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;
  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 0.0; 
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( a );
  ProteinNodes.push_back( new ProteinNode( Nodes.back() ) );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 1.0, 0.0; 
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( b );
  ProteinNodes.push_back( new ProteinNode( Nodes.back() ) );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 0.0; 
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( c );
  ProteinNodes.push_back( new ProteinNode( Nodes.back() ) );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 1.0; 
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( d );
  ProteinNodes.push_back( new ProteinNode( Nodes.back() ) );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0,-1.0; 
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( e );
  ProteinNodes.push_back( new ProteinNode( Nodes.back() ) );
 

  // Create materials
  srand(time(NULL));
  double epsilon = double(rand())/double(RAND_MAX);
  double sigma = double(rand())/double(RAND_MAX);
  
  ProteinLennardJones Mat1(sigma, epsilon);

  // Create potential bodies
  double SearchR = 3.0; 
  ProteinBody TestBody1(ProteinNodes, &Mat1, SearchR);

  // Check consistency
  TestBody1.checkConsistency(true); // verbose set to false, for more details set it to true

  // Print Body Energy
  TestBody1.compute(true, false, false);
  cout << "Protein body energy = " << TestBody1.energy() << endl;
  	
  return 0;
}



