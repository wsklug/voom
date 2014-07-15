#include "SCElastic.h"
#include "FVK.h"
#include <iostream>
#include "TriangleQuadrature.h"
#include "LoopShell.h"
#include "Node.h"
#include <vector>
#include "voom.h"
#include <tvmet/Vector.h>

#include <stdlib.h>
#include <time.h>
#include <stdio.h>

using namespace std;
using namespace voom;



int main()
{
  const bool checkNodalInfo = false;
  const bool checkMaterial  = true;
	
  // typedefs
  typedef vector< DeformationNode<3>* >	NodeContainer;
  typedef LoopShellShape::CornerValences  CornerValences;
	

  // data
  NodeContainer nodes;
//   typedef SCElastic Material_t;
//   Material_t bilayer(1.0, 1.0, 0.0);
  typedef FVK Material_t;
  FVK bilayer(1.0, 0.0, 0.0, 1.0, 0.4);

#define ICOS
#ifdef ICOS
  const int numberOfNodes = 9;
  const CornerValences V(5, 5, 5);
  char filename[] = "icos.info";
#else
  const int numberOfNodes = 12;
  const CornerValences V(6, 6, 6);
  char filename[] = "node.info";
#endif

  // open ifs to read nodal info.
  ifstream ifs;
  ifs.open(filename, ios::in);
  //ifs.open("icos.info", ios::in);
  if( !ifs.good() ) {
    cout << "can not open input file." << endl;
    exit(0);
  }
  for (int i = 0; i < numberOfNodes; i++) {
    unsigned id = 0;
    DeformationNode<3>::Point X;
    DeformationNode<3>::Point x;
    ifs >> id
	>> X(0) >> X(1) >> X(2) // read in ref positions
	>> x(0) >> x(1) >> x(2); // read in displacements
    cout << id << endl
	 << X(0) << " " <<X(1) << " " << X(2) << endl
	 << x(0) << " " <<x(1) << " " << x(2) << endl;
    x += X;

    NodeBase::DofIndexMap idx(3);
    idx[0] = 3*i; idx[1] = 3*i+1; idx[2] = 3*i+2;

    DeformationNode<3> * nd = new DeformationNode<3>(id,idx,X,x);

    nd -> setId(id);
    cout << "Read "<<i<<"th node from input: " << endl
	 << nd->id() << endl
	 << nd->position() << endl
	 << nd->point() << endl;
    nodes.push_back( nd );
  }

#ifndef ICOS
  for (NodeContainer::iterator i=nodes.begin(); i!=nodes.end(); i++){
    DeformationNode<3>::Point x = (*i)->position();
    x(2) = 1.0e-1*sqrt( x(0)*x(0) + x(1)*x(1) );
//     (*i)->setPosition( x );
    (*i)->setPoint( x );
  }
#endif

  srand(0);
  for (NodeContainer::iterator n=nodes.begin(); n!=nodes.end(); n++){
    for(int i=0; i<(*n)->dof(); i++) {
      (*n)->addPoint( i, 0.02*(0.5 - rand()/(double)(RAND_MAX)) );
    }
  }
  // end of initialization
  int qOrder;
  std::cout << "Input quadrature order:";
  cin >> qOrder;
  std::cout << std::endl;

  double pressure = 0.0;
  NodeBase::DofIndexMap idx(1); idx[0] = -1;
  MultiplierNode pressureNode(nodes.size(), idx);
  pressureNode.setPoint(pressure);
  double tension = 0.0;
  MultiplierNode tensionNode(nodes.size(), idx);
  tensionNode.setPoint(tension);
  double penaltyCoefficient = 0.0e4;

  // create elements
  TriangleQuadrature quad(qOrder);

  LoopShell<Material_t> shellElem(quad, bilayer, nodes, V,
				  &pressureNode, &tensionNode);
  cout << "Built shell element with " << numberOfNodes << " nodes." << endl;

  shellElem.checkConsistency();
  //	shellElem.checkPositoins();

//   int rank, dimension;
//   shellElem.rank(rank, dimension);
	
//   cout << "Rank of the element stiffness matrix is  "
//        << rank << " out of " << dimension << endl;;
	
	
  return 0;
}
