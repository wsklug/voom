#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include "Node.h"
#include "SemiflexibleGel.h"
#include "BrownianDynamics.h"
#include "BrownianRod.h"
#include "BrownianDynamics.h"
#include "Crosslink.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;



// To-do list:
// 
// 1. Write BrownianRod class to implement computation of Brownian
// kicks and viscous drags/mobilities.  Probably derive from Spring<N>.
//
// 2. Update Spring<N> to store tangent vector.  Make accessible for
// Angle and BrownianRod.
// 
// 3. Create Clock class to keep track of time.  Interface
// BrownianDyanmics with the Clock.
//
// 4. Derive class SinusoidalForce<N> from ConcentratedForce<N>, using
// Clock for time
// 
// 5. Implement crosslink class
//
// 6. Implement periodic boundary conditions
//
// 7. Create 1-D material class to generalize 1-D Hookean spring to
// nonlinear constitutive laws.
//

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
  // Parameters
  ////////////////////////////////////////////////////////////////////

  // Estimate kBond = E*A/L; E~1GPa, A~10nm^2, L~1000nm
  double kBond = 10.0; // pN/nm 

  //Estimate kAngle = E*I/L = 2kT\xi_p/L; kT=4.1pN-nm, \xi_p=10^4nm, L~100nm
  //make kAngle 100 times larger to simulate the rotation diffusion
  double kAngle = 820.0*1000; // pN-nm

  // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
  // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
  double viscosity = 1.0e-9; // MPa-s

  // drag = 2\pi*\eta_s*L; L=100nm, 
  // drag = 6*10^-7 pN-s/nm
  double drag = 6.0e-7; // pN-s/nm

  double kT = 4.1e0; // pN-nm

  double dt = 1.0e-7; // s

  double L = 2.0e3; // 1 micron = 10^3 nm

  //double angle = 30.0 * M_PI /180.0; // tilted 30 degree
  
  // create body 
  SemiflexibleGel<2> * gel = new SemiflexibleGel<2>;
  

  SemiflexibleGel<2>::DefNodeContainer nodes;

  int nRods = 2;
  int nNodesperRod = 5; //Number of nodes per rod
  int nCrosslink = 1; //Number of crosslinks
  int nNodes= nNodesperRod * nRods + nCrosslink * 2;
  for( int e=0; e<nRods; e++ ) {

    int a;
    unsigned id;
    BrownianNode<2>::Point X;
    NodeBase::DofIndexMap idx(2);
    SemiflexibleGel<2>::DefNode * nd;

    for ( int i=0; i<nNodesperRod; i++) {
      a=nNodesperRod*e+i;
      id = a;
      X = 0.0+L*i/(nNodesperRod-1), 0.0;
      idx[0]=2*a; idx[1]=2*a+1;
      nd = new SemiflexibleGel<2>::DefNode(id,idx,X,X);
      nd -> setId(id);
      nodes.push_back( nd );
    }
  }
  
  // define 2 nodes for Crosslink
  int  a = nNodes-2;
  unsigned id = a;
  BrownianNode<2>::Point X;
  NodeBase::DofIndexMap idx(2);
  SemiflexibleGel<2>::DefNode * nd;
  X = 0.0 , 0.0;
  idx[0]=2*a; idx[1]=2*a+1;
  nd = new SemiflexibleGel<2>::DefNode(id,idx,X,X);
  nd -> setId(id);
  nodes.push_back(nd);
  a = a+1;
  id = a;
  X = 0.0 , 0.0;
  idx[0]=2*a; idx[1]=2*a+1;
  nd = new SemiflexibleGel<2>::DefNode(id,idx,X,X);
  nd -> setId(id);
  nodes.push_back(nd);

  for( int e=0; e<nRods; e++ ) {
  
    // create a filament with springs and angle springs connecting the
    // series of nodes (here three nodes)
    SemiflexibleGel<2>::DefNodeContainer n(nNodesperRod);
    for( int i=0; i<nNodesperRod; i++) {
      n[i] = nodes[nNodesperRod*e+i];
    }
    gel->addFilament( n, kBond, kAngle, viscosity, kT, dt );

  }
  


  cout << "check point 1"<<endl;  
  // add a crosslink
  Crosslink<2> * crosslink = 
    new Crosslink<2> ( nodes[0], nodes[1], // 1A, 1B 
		       nodes[5], nodes[6], // 2A, 2B
		       nodes[10], nodes[11], // 1, 2
		       0.6, 0.4, 10.0);
    cout << "check point 2"<<endl;  
  gel->addCrosslink(crosslink);

  // create time integrator for Brownian dynamics
  cout << "check point 3"<<endl; 
  int printStride =  -1;
  int nSteps = 10000;
  int nPrintSteps = 100;

  BrownianDynamics * brownian = 
    new BrownianDynamics( nodes, printStride );

  brownian->pushBackBody( gel );

//   brownian->checkConsistency();

  ofstream pointTracker("txy.dat");
  ofstream forceTracker("fxy.dat");
  ofstream mobTracker("mxy.dat");

  
  for(int printStep=0; printStep<nPrintSteps; printStep++) {
    pointTracker << printStep*nSteps << ' ';
    for(int a=0; a<nNodes; a++) {
      pointTracker << nodes[a]->getPoint(0) << ' '
		   << nodes[a]->getPoint(1) << ' ';
      forceTracker << nodes[a]->getForce(0) << ' '
		   << nodes[a]->getForce(1) << ' ';
      mobTracker << nodes[a]->getMobility(0,0) << ' '
		 << nodes[a]->getMobility(1,1) << ' ';
    }
    pointTracker << std::endl;
    forceTracker << std::endl;
    mobTracker << std::endl;
   
    brownian->run( nSteps, dt );

    std::cout << "BrownianDynamics: step "<< printStep*nSteps
	      << std::setprecision( 16 ) 
	      << " | energy = "
	      << brownian->energy() << std::endl;
    char s[20];
    sprintf(s,"step%04d", printStep*nSteps);
    // print gel state to a file for visualization with ParaView
    std::cout << "Printing gel... ";
    std::cout.flush();
    gel->print(s);
    std::cout << "done. ";    
  }

  pointTracker.close();

  std::cout << "All done." << std::endl;
  return 0;

}
