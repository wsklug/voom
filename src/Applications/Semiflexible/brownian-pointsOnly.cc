#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include "Node.h"
#include "SemiflexibleGel.h"
#include "BrownianDynamics.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  if( argc != 1 ) {
      cout << "Usage: brownian" << endl;
      return(0);
  }

  // Alert user to OpenMP status
#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  // Initialize MPI if in use
  bool verbose=true;
#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif

  ////////////////////////////////////////////////////////////////////
  // Input section
  ////////////////////////////////////////////////////////////////////

  tvmet::Matrix<double, 3, 2> xn;
  xn = 
    -1.0, 1.0e-2,
    0.0, 0.0,
    1.0, 1.0e-2;
//   xn = 
//     0.0, 1.0,
//     0.0, 0.0,
//     1.0, 0.0;
  
  SemiflexibleGel<2>::DefNodeContainer nodes;

  int nNodes=100;
  for( int a=0; a<nNodes; a++ ) {

    unsigned id = a;
    DeformationNode<2>::Point X;
    DeformationNode<2>::Point x;
//     X = xn(a,0), xn(a,1);
    X = 0.0,0.0;

    NodeBase::DofIndexMap idx(2);
    idx[0] = 2*a; idx[1] = 2*a+1; 
    

    DeformationNode<2> * nd = new DeformationNode<2>(id,idx,X,X);
    
    nd -> setId(id);
    
    cout << "Read "<<a<<"th node from input: " << endl
	 << nd->id() << endl
	 << "X = " 
	 << std::setw(16) << nd->position()(0)
	 << std::setw(16) << nd->position()(1) 
	 << endl;
    nodes.push_back( nd );

  }

  // create body with one filament with three nodes
  SemiflexibleGel<2> * gel = new SemiflexibleGel<2>;

  // Estimate kBond = E*A/L; E~1GPa, A~10nm^2, L~100nm
  double kBond = 100.0; // pN/nm 

  //Estimate kAngle = E*I/L = 2kT\xi_p/L; kT=4.1pN-nm, \xi_p=10^4nm, L~100nm
  double kAngle = 820.0; // pN-nm

//   gel->addFilament( nodes, kBond, kAngle );

  // create time integrator for Brownian dynamics

  // drag = 2\pi*\eta_s*L; \eta_s = 10^3 Pa-s, L=100nm
  double drag = 6.0e-1; // pN-s/nm
  double kT = 4.1e0; // pN-nm
  double dt = 1.0e-4; // s
  int printStride =  -1;
  int nSteps = 1;
  int nPrintSteps = 10000;

  std::vector<NodeBase * > baseNodes( nodes.begin(), nodes.end() );
  BrownianDynamics * brownian = 
    new BrownianDynamics( baseNodes, drag, kT, dt, printStride );

//   brownian->pushBackBody( gel );

//   brownian->checkConsistency();

  ofstream pointTracker("txy.dat");
  for(int printStep=0; printStep<nPrintSteps; printStep++) {
    pointTracker << printStep*nSteps << ' ';
    for(int a=0; a<nNodes; a++) {
      pointTracker << nodes[a]->getPoint(0) << ' '
		   << nodes[a]->getPoint(1) << ' ';
    }
    pointTracker << std::endl;
    brownian->run( nSteps );

    std::cout << "BrownianDynamics: step "<< printStep*nSteps
	      << std::setprecision( 16 ) 
	      << " | residual norm = "
	      << blitz::max( blitz::abs( brownian->force() ) )
	      << " | energy = "
	      << brownian->energy() << std::endl;
    char s[20];
    sprintf(s,"step%04d", printStep*nSteps);
//     gel->print(s);
  }

  pointTracker.close();

  std::cout << "All done." << std::endl;
  return 0;
}
