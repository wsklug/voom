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
#include "Dirichlet.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;



// To-do list:
//
// 2. Update Spring<N> to store tangent vector.  Make accessible for
// Angle and BrownianRod.
// 
// 3. Create Clock class to keep track of time.  Interface
// BrownianDyanmics with the Clock.
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

   int nNodesperRod = 20; //Number of nodes per rod
  int nNodes= nNodesperRod;

  double L = 10.0; // micron

  double dL = L/(nNodesperRod-1); // 0.1 micron = 100 nm


  double kT = 4.1e-3; // pN-micron

  //Estimate kAngle = E*I/L = kT\xi_p/dL; kT=4.1pN-nm, \xi_p=10^4nm, dL~100nm
  //make kAngle 100 times larger to simulate the rotation diffusion
  double xi_p = 17.0; // microns
  double kAngle = kT*xi_p/dL; // pN-micron

  // Estimate kBond = EA/L; EA~ 12 EI/d^2 = 12*\xi_p*kT/d^2; d~8nm;
  // EA ~ 2/3 * 10^4 pN
  double r=3.5e-3; // radius in microns
  double kBond = 4.0*kAngle/(r*r); // pN/micron

  // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
  // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
//   double viscosity = 1.0e-9; // MPa-s
  double viscosity = 1.0; // pN-ms/micron^2

  // drag = 2\pi*\eta_s*L; L=100nm, 
  // drag = 6*10^-7 pN-s/nm

//   double dt = 1.0e-7; // s
  double dt = 1.0e-5; // ms

  //double angle = 30.0 * M_PI /180.0; // tilted 30 degree
  
  // create body 
  SemiflexibleGel<2> * gel = new SemiflexibleGel<2>;
  

  SemiflexibleGel<2>::DefNodeContainer nodes;

    int a;
    unsigned id;
    BrownianNode<2>::Point X(0.0);
    NodeBase::DofIndexMap idx(2);
    SemiflexibleGel<2>::DefNode * nd;
 
    ranlib::Normal<double> rng( 0.0, sqrt(dL/xi_p) );
    rng.seed((unsigned int)time(0));
    double dtheta_avg = 0.0;
    double dtheta2_avg = 0.0;

    blitz::Array<double,1> theta(nNodesperRod);
    theta = 0.0;
    for ( int a=0; a<nNodesperRod; a++) {      
      double dtheta = rng.random(); 
      theta(a) = theta( std::max(0,a-1) ) + dtheta;
      dtheta_avg += dtheta;
      dtheta2_avg += sqr(dtheta);
      std::cout << "a = " << a 
		<< " dtheta = " << dtheta*180.0/M_PI 
		<< " theta = " << theta(a)*180.0/M_PI 
		<< std::endl;
    }
    theta(nNodesperRod-1) = 0.0;

    // make the rod straight
    //theta = 0.0;

    for ( int a=0; a<nNodesperRod; a++) {
      id = a;
      idx[0]=2*a; idx[1]=2*a+1;
      nd = new SemiflexibleGel<2>::DefNode(id,idx,X,X);
      nd -> setId(id);
      nodes.push_back( nd );

      BrownianNode<2>::Point dX(0.0);
      dX = dL*cos(theta(a)), dL*sin(theta(a));
      X = X + dX;
    }
    dtheta_avg /= nNodesperRod;
    dtheta2_avg /= nNodesperRod;
    
    std::cout << " dtheta_avg = " << dtheta_avg 
	      << " dtheta2_avg = " << dtheta2_avg << std::endl;
     
    // Rotate filament so that ends are on x-axis 

    X = nodes.back()->point() - nodes.front()->point();
    double C = X(0)/norm2(X);
    double S = X(1)/norm2(X);
    std::cout << " X = " << X << std::endl 
	      << " L = " << norm2(X) <<  " cosine = " << C <<  " sine = " << S << std::endl;
    for ( int a=0; a<nNodesperRod; a++) {
      const BrownianNode<2>::Point & x = nodes[a]->point();
      nodes[a]->setPosition( 0,  C*x(0)+S*x(1) );
      nodes[a]->setPosition( 1, -S*x(0)+C*x(1) );
      nodes[a]->setPoint( nodes[a]->position() );
    }

    X = nodes.back()->point() - nodes.front()->point();
    std::cout << " X = " << X << std::endl;

    // create a filament with springs and angle springs connecting the
    // series of nodes (here three nodes)
    gel->addFilament( nodes, kBond, kAngle, viscosity, kT, dt );

    // pin first node
    gel->addConstraint( new Dirichlet( nodes[0], 0 ) );
    gel->addConstraint( new Dirichlet( nodes[0], 1 ) );
    // pin vertical displacement of last node
    gel->addConstraint( new Dirichlet( nodes[nNodesperRod-1], 1 ) );

//     gel->addConstraint( new Dirichlet( nodes[0], 0 ) );
//     // pin vertical displacements of all nodes
//     for(int a=0; a<nNodesperRod; a++) {
//       gel->addConstraint( new Dirichlet( nodes[a], 1 ) );
//     }


  // create time integrator for Brownian dynamics

  int printStride =  -1;
  int nSteps = 100000;
  int nPrintSteps = 100;

  BrownianDynamics * brownian = 
    new BrownianDynamics( nodes, printStride );

  brownian->pushBackBody( gel );

//   brownian->checkConsistency();

  ofstream avgTracker("avg.dat");
  ofstream pointTracker("txy.dat");
  ofstream forceTracker("fxy.dat");
  ofstream mobTracker("mxy.dat");

  double Lavg=0.0;

  for(int printStep=0; printStep<=nPrintSteps; printStep++) {
    pointTracker << printStep*nSteps << ' ';
    avgTracker << printStep*nSteps << ' ';

    for(int f=0; f<gel->filaments().size(); f++) {

      double xavg = 0.0;
      double yavg = 0.0;      
      double tavg = 0.0;      
      for(int a=0; a<gel->filament(f)->nodes.size()-1; a++) {
	SemiflexibleGel<2>::DefNode * na = gel->filament(f)->nodes[a];
	SemiflexibleGel<2>::DefNode * nb = gel->filament(f)->nodes[a+1];
	xavg += 0.5*( na->getPoint(0) + nb->getPoint(0) );
	yavg += 0.5*( na->getPoint(1) + nb->getPoint(1) );
	double dx = na->getPoint(0) - nb->getPoint(0);
	double dy = na->getPoint(1) - nb->getPoint(1);
	tavg +=  std::atan( dy/dx );
      }
      xavg /= gel->filament(f)->nodes.size()-1;
      yavg /= gel->filament(f)->nodes.size()-1;
      tavg /= gel->filament(f)->nodes.size()-1;

      avgTracker << xavg << ' '
		 << yavg << ' '
		 << tavg << ' ';

      SemiflexibleGel<2>::DefNode * first = gel->filament(f)->nodes.front();
      SemiflexibleGel<2>::DefNode * last  = gel->filament(f)->nodes.back();

      std::cout << " L = " << last->getPoint(0) - first->getPoint(0) << std::endl;

      Lavg += last->getPoint(0) - first->getPoint(0);

      std::cout << " Lavg = " << Lavg/(printStep+1) << std::endl;

      pointTracker << first->getPoint(0) << ' '
		   << first->getPoint(1) << ' ';
      forceTracker << first->getForce(0) << ' '
		   << first->getForce(1) << ' ';
      mobTracker << first->getDrag(0,0) << ' '
		 << first->getDrag(0,1) << ' '
		 << first->getDrag(1,0) << ' '
		 << first->getDrag(1,1) << ' ';

      pointTracker << last->getPoint(0) << ' '
		   << last->getPoint(1) << ' ';
      forceTracker << last->getForce(0) << ' '
		   << last->getForce(1) << ' ';
      mobTracker << last->getDrag(0,0) << ' '
		 << last->getDrag(0,1) << ' '
		 << last->getDrag(1,0) << ' '
		 << last->getDrag(1,1) << ' ';
    }
    pointTracker << std::endl;
    forceTracker << std::endl;
    mobTracker << std::endl;
    avgTracker << std::endl;
   
    std::cout << "BrownianDynamics: step "<< printStep*nSteps
	      << std::setprecision( 16 ) 
	      << " | energy = "
	      << brownian->energy() << std::endl;
    char s[20];
    sprintf(s,"step%04d", printStep*nSteps);
    // print gel state to a file for visualization with ParaView
    gel->print(s);

    brownian->run( nSteps, dt );

  }

  Lavg /= nPrintSteps+1;

  std::cout << "Lavg = " << Lavg << std::endl
	    << "dLavg = " << Lavg-L << std::endl;
  
  double dLavg = 1.645*0.5*L*L/(M_PI*M_PI*xi_p);
  std::cout << "Inextinsible Lavg = " << L - dLavg
	    << std::endl
	    << "Inextinsible dLavg = " << -dLavg
	    << std::endl;
  pointTracker.close();

  std::cout << "All done." << std::endl;
  return 0;

}
