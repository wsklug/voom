
#include "QuadQuadrature.h"
#include "ShapeQ4.h"
#include "CardiacPotential.h"
#include "Shape.h"
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

  VoltageNode<2>::PositionVector x;
  VoltageNode<2>::Point v;
  NodeBase::DofIndexMap idx(1);
  CardiacPotential::VoltageNodeContainer nodes;
  nodes.clear();
  idx[0]=0;

  int id;
 
  id=0;
  x[0]=0.0; x[1]=0.0; v=-20.0;
  nodes.push_back(new VoltageNode<2>(id,idx,x,v));

  id=0;
  x[0]=4.0; x[1]=5.0; v=-60.0;
  nodes.push_back(new VoltageNode<2>(id,idx,x,v));

  id=0;
  x[0]=6.0; x[1]=11.0; v=-80.0;
  nodes.push_back(new VoltageNode<2>(id,idx,x,v));

  id=0;
  x[0]=2.0; x[1]=6.0; v=-10.0;
  nodes.push_back(new VoltageNode<2>(id,idx,x,v));

  for (int i=0;i<nodes.size();i++) {
      double xcoord=nodes[i]->position(0);
      double ycoord=nodes[i]->position(1);
      double volt=nodes[i]->getPoint(0);
      cout << xcoord << " " <<ycoord <<" "<<volt<<endl;

  }

  QuadQuadrature q(2);
  CardiacPotential cp(q,nodes);
  cp.compute(0.0,true,0.0,false);


  return 0;
}
