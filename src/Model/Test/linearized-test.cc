// -*- C++ -*-
#include <vector>
#include <tvmet/Vector.h>
#include "TriangleQuadrature.h"
#include "LinearizedElement2D.h"
#include "LinearizedBody2D.h"
#include "Hookean2D.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "Node.h"
#include "Model.h"
#include "DirectLinearSolver.h"
using namespace std;
using namespace voom;

int main()
{
  const int nNodes = 13;
  double x[][2] = {
    {0.0  , 0.0	 },	// 0
    {1.0  , 0.0	 },	// 1
    {0.452, 0.593},	// 2
    {1.0  , 1.0	 },	// 3
    {0.0  , 1.0	 },	// 4
    {0.5  , 0.0	 },	// 5
    {0.75 , 0.25 },	// 6
    {0.25 , 0.25 },	// 7
    {1.0  , 0.5	 },	// 8
    {0.75 , 0.75 },	// 9
    {0.5  , 1.0	 },	// 10
    {0.25 , 0.75 },	// 11
    {0.0  , 0.5	 }	// 12
  };

  double u[][2] = {
    {0.00,	0.00},	// 0 
    {0.01,	0.03},	// 1 
    {0.00,	0.00},	// 2 
    {0.01,	0.05},	// 3 
    {0.00,	0.02},	// 4 
    {0.005,	0.015},	// 5 
    {0.00,	0.00},	// 6 
    {0.00,	0.00},	// 7 
    {0.01,	0.04},	// 8 
    {0.00,	0.00},	// 9 
    {0.005,	0.035},	// 10
    {0.00,	0.00},	// 11
    {0.00,	0.01}	// 12
  };

  bool constrained[][2] = {
    {true ,true },	// 0 
    {true ,true },	// 1 
    {false,false},	// 2 
    {true ,true },	// 3 
    {true ,true },	// 4 
    {true ,true },	// 5 
    {false,false},	// 6 
    {false,false},	// 7 
    {true ,true },	// 8 
    {false,false},	// 9 
    {true ,true },	// 10
    {false,false},	// 11
    {true ,true }	// 12
  };

  const int nElements = 4;
  const int nodesPerElement = 6;
  const int indx[][nodesPerElement] = {
    {0,1,2,5,6,7},
    {1,3,2,8,9,6},
    {0,2,4,7,11,12},
    {2,3,4,9,10,11}
  };

  std::vector< tvmet::Vector<int,6> > connectivity(nElements);
  for(int e=0; e<nElements; e++)
    for(int i=0; i<nodesPerElement; i++)
      connectivity[e](i) = indx[e][i];

  typedef voom::LinearizedElement2D<  
    voom::DeformationNode<2>, 
    voom::TriangleQuadrature, 
    voom::Hookean2D, 
    voom::ShapeTri6          > TriElement;

  typedef voom::LinearizedBody2D<  
    voom::DeformationNode<2>, 
    voom::TriangleQuadrature, 
    voom::Hookean2D, 
    voom::ShapeTri6          > TriBody;

  const int quadOrder = 2;

  bool verbose=false;
  
  Body::NodeContainer nodes;
  
  for(int a=0; a<nNodes; a++) {
    DeformationNode<2>::PositionVector xa;
    DeformationNode<2>::Point ua;
    for(int i=0; i<2; i++) {
      xa(i) = x[a][i];
      ua(i) = u[a][i];
    }
    NodeBase::DofIndexMap idx(2);
    idx[0] = 2*a; idx[1] = 2*a+1;
    DeformationNode<2> * n = new DeformationNode<2>(a, idx, xa, ua );
    nodes.push_back( n );
  }
  
  TriangleQuadrature quad(quadOrder);
  Hookean2D mat(Hookean2D::planeStrain, 10.0, 0.3);

  Body::ElementContainer elements;

  for(int e=0; e<nElements; e++) {
    TriElement::NodeContainer nodes_e(nodesPerElement);
    std::cout << "Constructing element "<<e<<" from nodes ";
    for(int n=0; n<nodesPerElement; n++) {
      nodes_e[n] = static_cast<DeformationNode<2>*>(nodes[connectivity[e][n]]);
      std::cout << nodes_e[n]->id() << " ";
    }
    std::cout << std::endl;
    TriElement * elem = new TriElement(quad,mat,nodes_e);
    elements.push_back( elem );
  }

  int nDof=2*nodes.size();
  for(int n=0; n<nNodes; n++) {
    for(int i=0; i<2; i++) {
      if( constrained[n][i] ) {
	MultiplierNodeConstraint<2>::Vector dir(0.0);
	dir(i) = 1.0;
	NodeBase::DofIndexMap idx(1);
	idx[0] = nDof;
	MultiplierNode * mnode = new MultiplierNode(nodes.size()-1, idx, 0.0);
	nodes.push_back(mnode);
	nDof++;

	MultiplierNodeConstraint<2> * mconstraint = 
	  new MultiplierNodeConstraint<2>(
	    static_cast<DeformationNode<2>*>(nodes[n]), mnode, dir, u[n][i]
	  );
	elements.push_back( mconstraint );
      }
    }
  }	   

  TriBody body(0, nodes, elements, TriBody::paraview);
  Model::BodyContainer bodies;
  bodies.push_back(&body);
  Model model( bodies );

  int ni=0;
  for(Body::ConstNodeIterator n=body.nodes().begin(); n!=body.nodes().end(); n++) {
    std::cout << "Node " << ni++ << ":" <<std::endl;
    for(NodeBase::DofIndexMap::const_iterator d=(*n)->index().begin();
	d != (*n)->index().end(); d++ ) {
      std::cout << std::setw(8) << (*d);
    }
    std::cout << std::endl;

  }

  int ei=0;
  for(Body::ConstElementIterator e=body.elements().begin(); e!=body.elements().end(); e++) {
    std::cout << "Element " << ei++ << ":" <<std::endl;
    for(ElementBase::DofIndexMap::const_iterator d=(*e)->index().begin();
	d != (*e)->index().end(); d++ ) {
      std::cout << std::setw(8) << (*d);
    }
    std::cout << std::endl;
  }

//   Storage solver;
//   solver.resize(model.dof());

//   model.computeAndAssemble(solver,true,true,true);

//   if(verbose) {
//     std::cout << "solver._E = "<<solver._E << std::endl
// 	    << "solver._DE = "<<solver._DE << std::endl
// 	    << "solver._DDE = "<<solver._DDE << std::endl;
//   }

//   if( model.checkConsistency() )
//     model.checkRank(model.dof());

  DirectLinearSolver linearSolver;
  linearSolver.resize(model.dof());
  linearSolver.solve(&model);
  model.computeAndAssemble(linearSolver,true,true,true);
  
  ei=0;
  for( Body::ConstElementIterator e=elements.begin(); e!=elements.end(); e++,ei++){
    TriElement* elem = dynamic_cast<TriElement*>(*e);
    if( !elem ) continue;

    std::cout << "Element " << ei << std::endl;
    const TriElement::QuadPointContainer & quadpts = elem->quadraturePoints();
    int pi=0;
    for( TriElement::ConstQuadPointIterator p=quadpts.begin(); 
	 p!=quadpts.end(); p++, pi++ ) {
      std::cout << "Point " << pi << std::endl;
      std::cout.precision(15);
      std::cout.setf(ios_base::fixed, ios_base::floatfield);
      std::cout <<"Strain = "<< p->material.strain() << std::endl
		<<"Stress = "<< p->material.stress() << std::endl;
    }
    std::cout << std::endl;
  }
  
  
  std::cout.precision(8);
  std::cout.setf(ios_base::fmtflags(0), ios_base::floatfield);

  std::cout <<setw(8) << "node" 
	    <<setw(34) << "x" 
	    <<setw(36) << "u"
	    <<setw(36) << "f"<<std::endl; 
  ni=0;
  for(Body::ConstNodeIterator n=nodes.begin(); ni<nNodes; n++,ni++) {
    DeformationNode<2> * nd = (DeformationNode<2>*)(*n);
    DeformationNode<2>::PositionVector x,u;
    x = nd->position();
    u = nd->point();
    std::cout<<std::setw(8)<<ni
	     <<std::setw(16)<<x(0)<<','
	     <<std::setw(16)<<x(1)<<std::setw(5)<<' '
	     <<std::setw(16)<<u(0)<<','
	     <<std::setw(16)<<u(1)
	     <<std::setw(16)<<linearSolver.gradient(2*ni)<<','
	     <<std::setw(16)<<linearSolver.gradient(2*ni+1)
	     <<std::endl;
  }

  std::cout << "Printing results to paraview file...";
  cout.flush();
  model.print("patchTest");
  std::cout << " done." << std::endl;
  return 0;

}
