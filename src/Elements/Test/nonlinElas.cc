// Test written by Melissa M. Gibbons
// for NonlinearElastic element class

#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
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

  // make 4 nodes for a tetrahedron
 
  std::vector< DeformationNode<3>* > defNodes;
  defNodes.reserve(4);

  typedef vector< DeformationNode<3>* >	NodeContainer;
  NodeContainer nodes;
 
  
  int dof=0;

  DeformationNode<3>::Point x;
  int id=1;
  x = 0.0, 0.0, 0.0;
  //x = -94.1547, 89.2191, -36.2134;
  //x/=10.0;
  NodeBase::DofIndexMap idx(3);
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
  nodes.push_back( n );
  defNodes.push_back( n );

  id=2;
  x = 2.0, 0.0, 0.0;
  //x = -94.1547, 91.4758, -50.6987;
  //x/=10.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* o = new DeformationNode<3>(id,idx,x);
  nodes.push_back( o );
  defNodes.push_back( o );

  id=3;
  x = 0.0, 2.0, 0.0;
  //x = -50.6987, 64.0328, -36.2134;
  //x/=10.0;
 
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* p = new DeformationNode<3>(id,idx,x);
  nodes.push_back( p );
  defNodes.push_back( p );

  id=4;
  x = 0.0, 0.0, 2.0;
  //x = -83.1061, 72.1094, -21.4441;
  //x/=10.0;
  
  for(int j=0; j<3; j++) idx[j]=dof++;
  DeformationNode<3>* q = new DeformationNode<3>(id,idx,x);
  nodes.push_back( q );
  defNodes.push_back( q );
  cout << "Number of nodes: " <<nodes.size() << endl;

  // connectivity for the tetrahedron
  vector< tvmet::Vector<int,4> > connectivities;
  tvmet::Vector<int, 4> c;
  connectivities.reserve(1);
  c[0]=1;
  c[0]=2;
  c[0]=3;
  c[0]=4;
  connectivities.push_back(c);

 //  typedef NonlinearElastic<
//     TetQuadrature,
//     StVenant,
//     ShapeTet4           > Tet4Element;

  // create material
  double rho = 500.0;
  double E = 200.0;//1.0e3;
  double nu = 0.4;
  typedef StVenant MaterialType;
  MaterialType protein( rho, E, nu );
  
  // create element
  TetQuadrature quad(1);

  typedef NonlinearElastic<TetQuadrature,StVenant,ShapeTet4> 
    ElementType;
  ElementType neElement(quad,protein,nodes);

  cout << "I initialized the element ok" << endl;

  //redundant
  //neElement.compute(true,true,false);
  
  cout << "Element energy = " << neElement.energy() << endl;
  cout << "Element deformation gradient (stored in material): " << neElement.quadraturePoints()[0].material.deformationGradient() << endl;
  neElement.checkConsistency();

  return 0;
}
