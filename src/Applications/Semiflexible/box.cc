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
#include <random/uniform.h>
#include <random/normal.h>
//#include "Dirichlet.h"

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

bool crossornot(BrownianNode<2>::Point x1a,
		BrownianNode<2>::Point x1b,
		BrownianNode<2>::Point x2a,
		BrownianNode<2>::Point x2b){
  double k1, k2, x, y;
  if (fabs(x1a(0)-x1b(0)) < 1.0e-6 && fabs(x2a(0)-x2b(0)) > 1.0e-6 ){
    k2 = (x2b(1)-x2a(1))/(x2b(0)-x2a(0));
    x = x1a(0);
    y = k2 *(x-x2b(0))+x2b(1);
    if (x1a(1)>x1b(1)) return (y < x1a(1) && y > x1b(1));
    else return (y > x1a(1) && y < x1b(1));
  }
  if (fabs(x1a(0)-x1b(0)) > 1.0e-6 && fabs(x2a(0)-x2b(0)) < 1.0e-6){
    k1 = (x1b(1)-x1a(1))/(x1b(0)-x1a(0));
    x = x2a(0);
    y = k1 *(x-x1b(0))+x1b(1);
    if ( x2a(1) > x2b(1) ) return ( y < x2a(1) && y > x2b(1));
    else return (y > x2a(1) && y < x2b(1));
  }
  if (fabs(x1a(0)-x1b(0)) < 1.0e-6 && fabs(x2a(0)-x2b(0)) < 1.0e-6){
    return false;
  }
  
  //insert main here
  double b1, b2;
  k1 = (x1b(1)-x1a(1))/(x1b(0)-x1a(0));
  k2 = (x2b(1)-x2a(1))/(x2b(0)-x2a(0));
  b1 = k1*x1a(0)-x1a(1);
  b2 = k2*x2a(0)-x2a(1);
  if (fabs(k1-k2) < 1.0e-6) return false; // two lines are parallel
  else {
    x = (b1-b2)/(k1-k2);
	y = k1*x - b1;
	if (min(x1a(0), x1b(0)) < x && x < max(x1a(0), x1b(0)) &&
        min(x2a(0), x2b(0)) < x && x < max(x2a(0), x2b(0)) &&
        min(x1a(1), x1b(1)) < y && y < max(x1a(1), x1b(1)) &&
        min(x2a(1), x2b(1)) < y && y < max(x2a(1), x2b(1)) )
		return true;
	else return false;
  }
}

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

  int nNodesperRod = 6; //Number of nodes per rod
  int nRods = 10;
  int nNodes = nNodesperRod * nRods;
  double Lbox = 10.0; 

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
  double kBond = (1.0e-4)*4.0*kAngle/(r*r); // pN/micron

  // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
  // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
//   double viscosity = 1.0e-9; // MPa-s
  double viscosity = 1.0; // pN-ms/micron^2

  // drag = 2\pi*\eta_s*L; L=100nm, 
  // drag = 6*10^-7 pN-s/nm

//   double dt = 1.0e-7; // s
  double dt = 1.0e-1; // ms

  //double angle = 30.0 * M_PI /180.0; // tilted 30 degree
  
  // create body 
  SemiflexibleGel<2> * gel = new SemiflexibleGel<2>;
  
  // all of the nodes from all of the filaments
  SemiflexibleGel<2>::DefNodeContainer nodes;


    
  for (int iRod=0; iRod<nRods; iRod++) {

    ranlib::Normal<double> rng( 0.0, sqrt(dL/xi_p) );
    ranlib::Uniform<double> rnu;
    rng.seed((unsigned int)time(0)+iRod);// has to seed every for loop otherwise the random numbers in different loops will be the same
    rnu.seed((unsigned int)time(0)+iRod);// has to seed every for loop otherwise the random numbers in different loops will be the same

    double dtheta_avg = 0.0;
    double dtheta2_avg = 0.0;
    
    blitz::Array<double,1> theta(nNodesperRod);


    theta(0) = rnu.random();
    theta(0) *= M_PI/2.0;
    std::cout << "a = " << 0 
		<< " dtheta = " << 0.0
		<< " theta = " << theta(0)*180.0/M_PI 
		<< std::endl;
    for ( int a=1; a<nNodesperRod; a++) {      
      double dtheta = rng.random(); 
      theta(a) = theta(a-1) + dtheta;
      dtheta_avg += dtheta;
      dtheta2_avg += sqr(dtheta);
      std::cout << "a = " << a 
		<< " dtheta = " << dtheta*180.0/M_PI 
		<< " theta = " << theta(a)*180.0/M_PI 
		<< std::endl;
    }
    //theta(nNodesperRod-1) = 0.0;
    
    // make the rod straight
    //theta = 0.0;


    // temporary container for the nodes of a particular filament
    SemiflexibleGel<2>::DefNodeContainer rod_nodes;

    // vector for nodal coordinates, initialized randomly for first node
    BrownianNode<2>::Point X(0.0);
    X = rnu.random()*Lbox, rnu.random()*Lbox;
    for ( int a=0; a<nNodesperRod; a++) {
      unsigned int id = a+nNodesperRod*iRod;
      NodeBase::DofIndexMap idx(2);
      idx[0]=2*id; idx[1]=2*id+1;
      SemiflexibleGel<2>::DefNode * nd 
	= new SemiflexibleGel<2>::DefNode(id,idx,X,X);
      
      // add to this filament's node list
      rod_nodes.push_back( nd );
      // add to global node list
      nodes.push_back( nd );
      
      BrownianNode<2>::Point dX(0.0);
      dX = dL*cos(theta(a)), dL*sin(theta(a));
      X = X + dX;
    }
    dtheta_avg /= nNodesperRod;
    dtheta2_avg /= nNodesperRod;
    
    std::cout << " dtheta_avg = " << dtheta_avg 
	      << " dtheta2_avg = " << dtheta2_avg << std::endl;
    
    // add the new filament to the gel body
    gel->addFilament( rod_nodes, kBond, kAngle, viscosity, kT, dt );


  }  // end of loop over filaments

  // Generate crosslinks
  //ofstream test1("test1.dat");
  SemiflexibleGel<2>::FilamentContainer filaments = gel->filaments();
  SemiflexibleGel<2>::FilamentIterator f1 = filaments.begin();
  for(; f1!=filaments.end(); f1++){
    //cout << "f1----" << std::endl;
    SemiflexibleGel<2>::FilamentIterator f2 = f1+1;
    for(; f2!=filaments.end(); f2++){
      //cout << "f2---" << std::endl;
      SemiflexibleGel<2>::RodContainer rods1 = (*f1)->rods;
      SemiflexibleGel<2>::RodContainer rods2 = (*f2)->rods;
      //cout << "**r1= " << rods1.size() << "**r2= " << rods2.size();
      SemiflexibleGel<2>::RodIterator r1 = rods1.begin();
      SemiflexibleGel<2>::RodIterator r2 = rods2.begin();
      for (; r1!=rods1.end(); r1++){
	//cout << "r1--" << std::endl;
	r2 = rods2.begin();
	for (; r2!=rods2.end();r2++){
	  //cout << "r2-" << std::endl;
	  SemiflexibleGel<2>::DefNodeContainer rodnodes1 = (*r1)->getNodes();
	  SemiflexibleGel<2>::DefNodeContainer rodnodes2 = (*r2)->getNodes();
	  const BrownianNode<2>::Point & X1A = rodnodes1[0]->position();
	  const BrownianNode<2>::Point & X1B = rodnodes1[1]->position();
	  const BrownianNode<2>::Point & X2A = rodnodes2[0]->position();
	  const BrownianNode<2>::Point & X2B = rodnodes2[1]->position();
	  /*SemiflexibleGel<2>::DefNodeContainer rodnodes1 = (*r1) -> getNodes();
	  X1A = rodnodes1[0]->point();
	  X1B = rodnodes1[1]->point();
	  SemiflexibleGel<2>::DefNodeContainer rodnodes2 = (*r2) -> getNodes();
	  X2A = rodnodes2[0]->point();
	  X2B = rodnodes2[1]->point();*/

	  /*
	  std::cout << "points tested " << " X1A= " << X1A(0)  
			<< " X1B= " << X1B(0) << std:: endl
                                          << " X2A= " << X2A(0)  
                                          << " X2B= " << X2B(0) << std::endl;  
	  std:cout << crossornot(X1A,X1B,X2A,X2B) << std::endl;*/
	  if(crossornot(X1A,X1B,X2A,X2B)){
	    Crosslink<2> * crosslink = new Crosslink<2> ( rodnodes1[0], rodnodes1[1], rodnodes2[0], rodnodes2[1], 0.6, 0.4, 10.0);
	    gel->addCrosslink(crosslink);
        std::cout << "crosslink created " << " X1A= " << X1A  
			<< " X1B= " << X1B << std:: endl
                                          << " X2A= " << X2A  
                                          << " X2B= " << X2B << std::endl;  
	  }
	}
      }
    }
  }
  // create time integrator for Brownian dynamics
  
  int printStride =  -1;
  int nSteps = 100;
  int nPrintSteps = 50;
  
  BrownianDynamics * brownian = 
    new BrownianDynamics( nodes, printStride );
  
  brownian->pushBackBody( gel );
  
  //   brownian->checkConsistency();

  for(int printStep=0; printStep<=nPrintSteps; printStep++) {
   
    std::cout << "BrownianDynamics: step "<< printStep*nSteps
	      << std::setprecision( 16 ) 
	      << " | energy = "
	      << brownian->energy() << std::endl;
    char s[20];
    sprintf(s,"step%04d", printStep*nSteps);
    // print gel state to a file for visualization with ParaView
    cout << "check point 0" << endl;
    gel->print(s);

    brownian->run( nSteps, dt );
    
  } // end print step loop

  std::cout << "All done." << std::endl;
  return 0;

}
