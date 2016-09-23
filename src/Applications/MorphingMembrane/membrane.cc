#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "FVK.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "C0MembraneBody.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "Lbfgsb.h"
#include "Model.h"
// #include "mesh.h"
//#include "mesh-order2.h"
// #include "mesh-order3.h"
#include "mesh-order4.h"

// #define NODENUMBER 42
// #define ELEMNUMBER 80

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{
	
  typedef C0MembraneBody<TriangleQuadrature,FVK,ShapeTri3> BodyType;
  //
  // create vector of nodes

  BodyType::NodeContainer nodes;
  for ( int a = 0; a < NumberOfNodes; a++){
    int id=-1;
    DeformationNode<3>::Point x;
    id = a;
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) {
      x(j) = CoordinateArray[a][j];
      idx[j] = 3*a+j;
    }
    nodes.push_back(new DeformationNode<3>(id,idx,x));
  }
  cout << nodes.size() << endl;
  for(int i=0; i<NumberOfNodes; i++) {
    cout << setw(12) << nodes[i]->id()
	 << setw(20) << nodes[i]->getPoint(0)
	 << setw(20) << nodes[i]->getPoint(1)
	 << setw(20) << nodes[i]->getPoint(2) 
	 << endl;
  }	
  //
  // create connectivities
  BodyType::ConnectivityContainer connectivities;
  BodyType::ElementConnectivity c(3);
  for (int e = 0; e < NumberOfElements; e++){
    for(int n=0; n<3; n++) {
      c[n] = Connectivity[e][n];
    }
    connectivities.push_back(c);
  }
  cout << connectivities.size() << endl;
  for(int i=0; i<NumberOfElements; i++) {
    cout << setw(12) << connectivities[i][0] 
	 << setw(12) << connectivities[i][1] 
	 << setw(12) << connectivities[i][2]
	 << endl;
  }
	
  const double kC=0.0;
  const double kG=0.0;
  const double C0=0.0;
  const double E=1.0;
  const double nu=0.3;
  FVK material(kC,kG,C0,E,nu);
  std::cout << "FVK Material has been created." << std::endl;
	
  // create Body
  double pressure=1.0;
  int quadOrder = 1;
  BodyType bd(material, connectivities, nodes, quadOrder, pressure);

  cout << "Created a body." << endl;
  bd.setOutput(Body::paraview);
  bd.compute(true,true, false);
  cout << "volume of current body = " << bd.volume() << endl;

  Model::BodyContainer bodies;
  bodies.push_back( &bd );
  Model model(bodies);
//   m.checkConsistency(true, false);
//   m.checkRank(m.dof()-6, true);

  // set up bounds for solver
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;
  for (int a=0; a<nodes.size(); a++) {  
    double x = nodes[a]->getPoint(0);
    double y = nodes[a]->getPoint(1);
    double eps=1.0e-6;
    if( std::abs(x) < eps || std::abs(x-1.0) < eps || std::abs(y) < eps || std::abs(y-1.0) < eps ) {
      std::cout << "Fixing node " << a << " at " 
		<< std::setw(24) << nodes[a]->getPoint(0)
		<< std::setw(24) << nodes[a]->getPoint(1)
		<< std::setw(24) << nodes[a]->getPoint(2) << std::endl;
      for(int i=0; i<3; i++) {
	nbd(3*a+i) = 2;
	hi(3*a+i) = lo(3*a+i) = nodes[a]->getPoint(i);
      }
    }
  }

  // Initialize nonlinear solver.
  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  double pentol=1.0e-4;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  std::cout << "Input iprint: " << iprint << std::endl
	    << "Input factr: " << factr << std::endl
	    << "Input pgtol: " << pgtol << std::endl
	    << "Input m: " << m << std::endl
	    << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);
  solver.setBounds(nbd,lo,hi);
  model.print("InitialState");
  solver.solve( &model );
  model.print("FinalState");

  return 0;
}

