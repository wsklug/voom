// Test for birck elements written by Luigi Perotti from
// Melissa M. Gibbons test for NonlinearElastic element class

#include "BrickQuadrature.h"
#include "ShapeBrick6.h"
#include "Shape.h"
#include "BrickElement.h"
#include "CompNeoHookean.h"
#include <vector>
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

using namespace std;
using namespace voom;

int main()
{ 
  // Make 6 nodes for a 6 node brick element
  std::vector< DeformationNode<3>* > Nodes;
  Nodes.reserve(6);

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
  x = 0.0, 1.0, 1.0; 
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( e );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 1.0; 
  DeformationNode<3>* f = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( f );

  

  // create material
  double rho = 500.0;
  double E = 200.0;
  double nu = 0.4;
  CompNeoHookean protein( rho, E, nu );
   


  // create element
  BrickQuadrature quad(1);
  ShapeBrick6 shape;
 
  BrickElement Brick(quad, protein, Nodes, shape);
  
  
  cout << "I initialized the element ok" << endl;
  
  cout << "Element energy = " << Brick.energy() << endl;
  cout << "Element deformation gradient, QP1: " << Brick.quadraturePoints()[0].material.deformationGradient() << endl;
  cout << "Element 1st PK stress, QP1: " <<  Brick.quadraturePoints()[0].material.piolaStress() << endl;

  Brick.checkConsistency();
  
  return 0;
}
