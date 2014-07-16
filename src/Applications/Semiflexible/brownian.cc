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
#include "Crosslink.h"
#include "Motor.h"
#include "GelOutput.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[]) {
  // Estimate kBond = E*A/L; E~1GPa, A~10nm^2, L~1000nm
  double kBond = 10.0; // pN/nm 

  //Estimate kAngle = E*I/L = 2kT\xi_p/L; kT=4.1pN-nm, \xi_p=10^4nm, L~100nmfor
  //make kAngle 100 times larger to simulate the rotation diffusion
  double kAngle = 820.0; // pN-nm

  // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
  // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
  double viscosity = 1.0e-9; // MPa-s

  // drag = 2\pi*\eta_s*L; L=100nm, 
  // drag = 6*10^-7 pN-s/nm
  double drag = 6.0e-7; // pN-s/nm

  double kT = 4.1e0; // pN-nm

  double dt = 1.0e-6; // s

  double L = 1.0e3; // 1 micron = 10^3 nm

  double kcl = -1.0; // pN/nm

  double kC = 5.0e4;

  tvmet::Vector<double,2> syssize;
  syssize =15.0,15.0;

  SemiflexibleGel<2>::DefNodeContainer dNodes;
  
  // create body 
//   SemiflexibleGel<2> * gel = new SemiflexibleGel<2>(dNodes,syssize,.25,9,.15,kAngle,viscosity,kT,dt,kcl,0.0,kC,3);
  SemiflexibleGel<2> * gel = new SemiflexibleGel<2>(dNodes,syssize,.75,10,.12,kBond,kAngle,viscosity,kT,dt,kcl,0.0);

  GelOutput<2> * gout = new GelOutput<2>();

  int nNodes = dNodes.size();
  
//     Motor<2> * motor = new Motor<2>(dNodes[0]->point(),.0075);
// 
//     gel->attachMotor(motor,(gel->filaments())[0],(gel->filaments())[2]);
// 
//     motCont.push_back(motor);
//   
//   //ugly hack to deal with static template problem //
//     motor->setFVParams(500.0,2.0,75.0,-200.0,200.0);
//     motor->setDetachParams(.1,2.0);
  // create time integrator for Brownian dynamics
  int printStride =  -1;
  int nSteps = 1000;
  int nPrintSteps = 1000;

  BrownianDynamics * brownian = 
    new BrownianDynamics( dNodes, printStride );

//   BrownianDynamics * brownian = 
//     new BrownianDynamics( dNodes, printStride );

  brownian->pushBackBody( gel );

//   brownian->checkConsistency();

//   ofstream pointTracker("txy.dat");
//   ofstream forceTracker("fxy.dat");
//   ofstream mobTracker("mxy.dat");
//   gel->print("step0");
  for(int printStep=0; printStep<nPrintSteps; printStep++) {
//     pointTracker << printStep*nSteps << ' ';
//     for(int a=0; a<nNodes; a++) {
//       pointTracker << dNodes[a]->getPoint(0) << ' '
// 		   << dNodes[a]->getPoint(1) << ' ';
//       forceTracker << dNodes[a]->getForce(0) << ' '
// 		   << dNodes[a]->getForce(1) << ' ';
//       mobTracker << dNodes[a]->getMobility(0,0) << ' '
// 		 << dNodes[a]->getMobility(1,1) << ' ';
//     }
/*    pointTracker << std::endl;
    forceTracker << std::endl;
    mobTracker << std::endl;*/
   
    brownian->run( nSteps, dt );
    std::cout << "t = " << dt*nSteps*(printStep+1) << std::endl;
//     std::cout << "BrownianDynamics: step "<< printStep*nSteps
// 	      << std::setprecision( 16 ) 
// 	      << " | energy = "
// 	      << brownian->energy() << std::endl;
    char s[20];
    sprintf(s,"step%04d", printStep*nSteps);
    // print gel state to a file for visualization with ParaView
//     std::cout << "Printing gel... ";
    std::cout.flush();
    (*gout)(gel,s);
/*    std::cout << "done. ";*/    
  }

//   pointTracker.close();

  std::cout << "All done." << std::endl;
  return 0;
}
