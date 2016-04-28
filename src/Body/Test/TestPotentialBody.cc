#include <string>
#include <iostream>
#include <vector>


#include "Node.h"
#include "SpringPotential.h"
#include "LennardJones.h"
#include "LennardJonesFT.h"
#include "Morse.h"

#include "PotentialBody.h"

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

 

  // Create materials
  srand(time(NULL));
  double epsilon = double(rand())/double(RAND_MAX);
  double sigma = double(rand())/double(RAND_MAX);
  double Rshift = double(rand())/double(RAND_MAX);
  LennardJones Mat1(sigma, epsilon);
  SpringPotential Mat2(sigma, epsilon);	
  LennardJonesFT Mat3(sigma, epsilon, Rshift);
  Morse Mat4(sigma, epsilon, Rshift);

  // Create potential bodies
  double SearchR = 3.0; 
  PotentialBody TestBody1(&Mat1, Nodes, SearchR);
  PotentialBody TestBody2(&Mat2, Nodes, SearchR);
  PotentialBody TestBody3(&Mat3, Nodes, SearchR);
  PotentialBody TestBody4(&Mat4, Nodes, SearchR);
  
  // Check consistency
  TestBody1.checkConsistency(false); // verbose set to false, for more details set it to true
  TestBody2.checkConsistency(false);
  TestBody3.checkConsistency(false);
  TestBody4.checkConsistency(false);
  
  // Print Body Energy
  TestBody1.compute(true, false, false);
  cout << "Potential body energy = " << TestBody1.totalStrainEnergy() << endl;
  TestBody2.compute(true, false, false);
  cout << "Potential body energy = " << TestBody2.totalStrainEnergy() << endl;
  TestBody3.compute(true, false, false);
  cout << "Potential body energy = " << TestBody3.totalStrainEnergy() << endl;
  TestBody4.compute(true, false, false);
  cout << "Potential body energy = " << TestBody4.totalStrainEnergy() << endl;
	
  return 0;
}



