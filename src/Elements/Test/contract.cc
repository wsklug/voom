// Test written by Mainak
// for FiniteContraction element class

#include "TetQuadrature.h"
//#include "ShapeTet10.h"
#include "ShapeTet4CP.h"
#include "StVenant.h"
#include "Shape.h"
#include "FiniteContraction.h"
#include "ContractionWrapper.h"
#include <vector>
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <time.h>
#include <iostream>
#include <math.h>

#include <cstdlib>

using namespace std;
using namespace voom;

int main()
{ 

  // make 4 nodes for a tetrahedron
 
  std::vector< DeformationNode<3>* > defNodes;
  defNodes.reserve(4);

  std::vector< ScalarFieldNode<3>* > voltNodes;
  voltNodes.reserve(4);

//  typedef vector< DeformationNode<3>* >	NodeContainer;
//  NodeContainer nodes;
 
  int dof=0;

//  std::vector< DeformationNode<3>::Point* > xVec;
//  xVec.reserve(4);
  double volts = 0.0;

  
  std::vector<Vector3D> X(4);
  X[0] = 0.0, 0.0, 0.0;
  X[1] = 2.0, 0.0, 0.0;
  X[2] = 0.0, 2.0, 0.0;
  X[3] = 0.0, 0.0, 2.0;

  srand(time(0));

  for(int id=0; id<4; id++) {

    Vector3D x; // current (deformed) position
    x = X[id];

    // perturb current positions by random displacement
    for(int i=0; i<3; i++) x(i) += 0.1*rand()/((double)RAND_MAX);

    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    defNodes.push_back( new DeformationNode<3>(id,idx,X[id],x) );
    voltNodes.push_back( new ScalarFieldNode<3>(id,idx,X[id],volts) );
    
  }
      
  cout << "Number of nodes: " <<defNodes.size() << endl;

  // connectivity for the tetrahedron
  vector< tvmet::Vector<int,4> > connectivities;
  tvmet::Vector<int, 4> c;
  connectivities.reserve(1);
  c[0] = 1;
    c[1] = 2;
    c[2] = 3;
    c[3] = 4;
//    c[4] = 5;
//    c[5] = 6;
//    c[6] = 7;
//    c[7] = 8;
//    c[8] = 9;
//    c[9] = 10;
  connectivities.push_back(c);

  // create material
  double rho = 500.0;
  double E = 200.0;//1.0e3;
  double nu = 0.4;
  typedef StVenant MaterialType;
  MaterialType passiveMuscle( rho, E, nu );
  
  typedef ContractionWrapper<StVenant> ContractionWrapMat;
  ContractionWrapMat activeMuscle(passiveMuscle);

  // create element
  TetQuadrature quad(2);

//  typedef NonlinearElastic<TetQuadrature,StVenant,ShapeHex8> 
//    ElementType;
//  ElementType neElement(quad,cardiacElement,nodes);

//  template<class DefQuadrature_t, class Material_t, class DefShape_t, 
//    		   class VoltQuadrature_t, class EpMaterial_t, class VoltShape_t>
  
  typedef FiniteContraction<TetQuadrature, ContractionWrapMat, ShapeTet4,
  TetQuadrature, ContractionWrapMat, ShapeTet4>
  ElementType;
  
  ElementType neElement(quad,activeMuscle,defNodes,quad,activeMuscle,voltNodes);
  
  cout << "I  initialized the element ok" << endl;

  //redundant
  //neElement.compute(true,true,false);
  
//  cout << "Element energy = " << neElement.energy() << endl;
//  cout << "Element deformation gradient, QP1: " << neElement.quadraturePoints()[0].material.deformationGradient() << endl;
//  cout << "Element deformation gradient, QP2: " << neElement.quadraturePoints()[1].material.deformationGradient() << endl;
//  cout << "Element deformation gradient, QP3: " << neElement.quadraturePoints()[2].material.deformationGradient() << endl;
//  cout << "Element deformation gradient, QP4: " << neElement.quadraturePoints()[3].material.deformationGradient() << endl;
//  cout << "Element 1st PK stress, QP1: " << neElement.quadraturePoints()[0].material.piolaStress() << endl;
//  cout << "Element 1st PK stress, QP2: " << neElement.quadraturePoints()[1].material.piolaStress() << endl;
//  cout << "Element 1st PK stress, QP3: " << neElement.quadraturePoints()[2].material.piolaStress() << endl;
//  cout << "Element 1st PK stress, QP4: " << neElement.quadraturePoints()[3].material.piolaStress() << endl;
//  cout << "Element Cauchy stress, QP1: " << neElement.quadraturePoints()[0].material.cauchyStress() << endl;
//  cout << "Element Cauchy stress, QP2: " << neElement.quadraturePoints()[1].material.cauchyStress() << endl;
//  cout << "Element Cauchy stress, QP3: " << neElement.quadraturePoints()[2].material.cauchyStress() << endl;
//  cout << "Element Cauchy stress, QP4: " << neElement.quadraturePoints()[3].material.cauchyStress() << endl;
  neElement.checkConsistency();

  return 0;
}
