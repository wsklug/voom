// Test written by Melissa M. Gibbons
// for NonlinearElastic element class
// tested as a quadratic tetrahedron

#include "TetQuadrature.h"
#include "ShapeTet10.h"
#include "StVenant.h"
#include "Shape.h"
#include "NonlinearElastic.h"
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

  // make 10 nodes for a tetrahedron
 
  std::vector< DeformationNode<3>* > defNodes;
  defNodes.reserve(10);

  typedef vector< DeformationNode<3>* >	NodeContainer;
  NodeContainer nodes;
 
  int dof=0;

  DeformationNode<3>::Point x;
  int id=1;
  x = 0.0, 0.0, 0.0;
 
  NodeBase::DofIndexMap idx(3);
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
  nodes.push_back( n );
  defNodes.push_back( n );

  id=2;
  x = 2.0, 0.0, 0.0;

  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* o = new DeformationNode<3>(id,idx,x);
  nodes.push_back( o );
  defNodes.push_back( o );

  id=3;
  x = 0.0, 2.0, 0.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* p = new DeformationNode<3>(id,idx,x);
  nodes.push_back( p );
  defNodes.push_back( p );

  id=4;
  x = 0.0, 0.0, 2.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* q = new DeformationNode<3>(id,idx,x);
  nodes.push_back( q );
  defNodes.push_back( q );

  id=5;
  x = 1.0, 0.0, 0.0;
  
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* r = new DeformationNode<3>(id,idx,x);
  nodes.push_back( r );
  defNodes.push_back( r );

  id=6;
  x = 1.0, 1.0, 0.0;
  
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* s = new DeformationNode<3>(id,idx,x);
  nodes.push_back( s );
  defNodes.push_back( s );

  id=7;
  x = 0.0, 1.0, 0.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* t = new DeformationNode<3>(id,idx,x);
  nodes.push_back( t );
  defNodes.push_back( t );

  id=8;
  x = 0.0, 0.0, 1.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* u = new DeformationNode<3>(id,idx,x);
  nodes.push_back( u );
  defNodes.push_back( u );

  id=9;
  x = 1.0, 0.0, 1.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* v = new DeformationNode<3>(id,idx,x);
  nodes.push_back( v );
  defNodes.push_back( v );

  id=10;
  x = 0.0, 1.0, 1.0;
  
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* w = new DeformationNode<3>(id,idx,x);
  nodes.push_back( w );
  defNodes.push_back( w );

  cout << "Number of nodes: " <<nodes.size() << endl;

  // connectivity for the tetrahedron
  vector< tvmet::Vector<int,10> > connectivities;
  tvmet::Vector<int, 10> c;
  connectivities.reserve(1);
  c[0] = 1;
  c[1] = 2;
  c[2] = 3;
  c[3] = 4;
  c[4] = 5;
  c[5] = 6;
  c[6] = 7;
  c[7] = 8;
  c[8] = 9;
  c[9] = 10;
  connectivities.push_back(c);

  // create material
  double rho = 500.0;
  double E = 200.0;//1.0e3;
  double nu = 0.4;
  typedef StVenant MaterialType;
  MaterialType protein( rho, E, nu );
  
  // create element
  TetQuadrature quad(2);

  typedef NonlinearElastic<TetQuadrature,StVenant,ShapeTet10> 
    ElementType;
  ElementType neElement(quad,protein,nodes);

  cout << "I initialized the element ok" << endl;

  //redundant
  //neElement.compute(true,true,false);
  
  cout << "Element energy = " << neElement.energy() << endl;
  cout << "Element deformation gradient, QP1: " << neElement.quadraturePoints()[0].material.deformationGradient() << endl;
  cout << "Element deformation gradient, QP2: " << neElement.quadraturePoints()[1].material.deformationGradient() << endl;
  cout << "Element deformation gradient, QP3: " << neElement.quadraturePoints()[2].material.deformationGradient() << endl;
  cout << "Element deformation gradient, QP4: " << neElement.quadraturePoints()[3].material.deformationGradient() << endl;
  cout << "Element 1st PK stress, QP1: " << neElement.quadraturePoints()[0].material.piolaStress() << endl;
  cout << "Element 1st PK stress, QP2: " << neElement.quadraturePoints()[1].material.piolaStress() << endl;
  cout << "Element 1st PK stress, QP3: " << neElement.quadraturePoints()[2].material.piolaStress() << endl;
  cout << "Element 1st PK stress, QP4: " << neElement.quadraturePoints()[3].material.piolaStress() << endl;
  cout << "Element Cauchy stress, QP1: " << neElement.quadraturePoints()[0].material.cauchyStress() << endl;
  cout << "Element Cauchy stress, QP2: " << neElement.quadraturePoints()[1].material.cauchyStress() << endl;
  cout << "Element Cauchy stress, QP3: " << neElement.quadraturePoints()[2].material.cauchyStress() << endl;
  cout << "Element Cauchy stress, QP4: " << neElement.quadraturePoints()[3].material.cauchyStress() << endl;
  neElement.checkConsistency();

  return 0;
}
