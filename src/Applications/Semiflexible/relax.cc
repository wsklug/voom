#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include "Node.h"
#include "SemiflexibleGel.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "BrownianRod.h"
#include <random/uniform.h>
#include <random/normal.h>
//#include "Dirichlet.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "GelOutput.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;


// bool crossornot(BrownianNode<2>::Point x1a,
// 		BrownianNode<2>::Point x1b,
// 		BrownianNode<2>::Point x2a,
// 		BrownianNode<2>::Point x2b){
//   double k1, k2, x, y;
//   if (fabs(x1a(0)-x1b(0)) < 1.0e-6 && fabs(x2a(0)-x2b(0)) > 1.0e-6 ){
//     k2 = (x2b(1)-x2a(1))/(x2b(0)-x2a(0));
//     x = x1a(0);
//     y = k2 *(x-x2b(0))+x2b(1);
//     if (x1a(1)>x1b(1)) return (y < x1a(1) && y > x1b(1));
//     else return (y > x1a(1) && y < x1b(1));
//   }
//   if (fabs(x1a(0)-x1b(0)) > 1.0e-6 && fabs(x2a(0)-x2b(0)) < 1.0e-6){
//     k1 = (x1b(1)-x1a(1))/(x1b(0)-x1a(0));
//     x = x2a(0);
//     y = k1 *(x-x1b(0))+x1b(1);
//     if ( x2a(1) > x2b(1) ) return ( y < x2a(1) && y > x2b(1));
//     else return (y > x2a(1) && y < x2b(1));
//   }
//   if (fabs(x1a(0)-x1b(0)) < 1.0e-6 && fabs(x2a(0)-x2b(0)) < 1.0e-6){
//     return false;
//   }
//   
//   double b1, b2;
//   k1 = (x1b(1)-x1a(1))/(x1b(0)-x1a(0));
//   k2 = (x2b(1)-x2a(1))/(x2b(0)-x2a(0));
//   b1 = k1*x1a(0)-x1a(1);
//   b2 = k2*x2a(0)-x2a(1);
//   if (fabs(k1-k2) < 1.0e-6) return false; // two lines are parallel
//   else {
//     x = (b1-b2)/(k1-k2);
// 	y = k1*x - b1;
// 	if (min(x1a(0), x1b(0)) < x && x < max(x1a(0), x1b(0)) &&
//         min(x2a(0), x2b(0)) < x && x < max(x2a(0), x2b(0)) &&
//         min(x1a(1), x1b(1)) < y && y < max(x1a(1), x1b(1)) &&
//         min(x2a(1), x2b(1)) < y && y < max(x2a(1), x2b(1)) )
// 		return true;
// 	else return false;
//   }
// }

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

  int nNodesperRod = 7; //Number of nodes per rod
  int nRods = 10;
  int nNodes = nNodesperRod * nRods;
  double Lbox = 8.0; 
  tvmet::Vector<double,2> syssize;
  syssize = Lbox,Lbox;

  double L = 2.0; // micron

  double dL = L/(nNodesperRod-1); // 0.1 micron = 100 nm


  double kT = 4.1e-3; // pN-micron

  int fitOrder = 3;

  //Estimate kAngle = E*I/L = kT\xi_p/dL; kT=4.1pN-nm, \xi_p=10^4nm, dL~100nm
  //make kAngle 100 times larger to simulate the rotation diffusion
  double xi_p = 17.0; // microns
  double kAngle = kT*xi_p/dL; // pN-micron

  double L_p = 15.0; // persistance length, microns
  double kC = L_p * kT; 

  // Estimate kBond = EA/L; EA~ 12 EI/d^2 = 12*\xi_p*kT/d^2; d~8nm;
  // EA ~ 2/3 * 10^4 pN
  double r=3.5e-3; // radius in microns
  double kBond = (1.0e-4)*4.0*kAngle/(r*r); // pN/micron

  // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
  // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
//   double viscosity = 1.0e-9; // MPa-s
  double viscosity = 0.0;//1.0; // pN-ms/micron^2

  // drag = 2\pi*\eta_s*L; L=100nm, 
  // drag = 6*10^-7 pN-s/nm

//   double dt = 1.0e-7; // s
  double dt = 2.5e-5; // ms

  double kcl = 1.0;
  double shrStart = 0.0;
  double shrEnd = 5.0;
  double shrStep = .5;

  double filDens = 15.0;

  char shrfilename[50];
  sprintf(shrfilename,"energyvsshear.dat");

  for(double shr = shrStart; shr <= shrEnd; shr += shrStep) {

  //double angle = 30.0 * M_PI /180.0; // tilted 30 degree
    SemiflexibleGel<2>::DefNodeContainer nodes;
  // create body 
    SemiflexibleGel<2> * gel = new SemiflexibleGel<2>(nodes,syssize,filDens,nNodesperRod,dL,kBond,kAngle,viscosity,kT,dt,kcl,shr);
  
  // all of the nodes from all of the filaments



  //////////////////////////////////////////////////////////////////////////
  //
  // Build a network by randomly placing rods with fluctuations built in
  //
  //////////////////////////////////////////////////////////////////////////

 //  for (int iRod=0; iRod<nRods; iRod++) {
// 
//     ranlib::Normal<double> rng( 0.0, sqrt(dL/xi_p) );
//     ranlib::Uniform<double> rnu;
//     rng.seed((unsigned int)time(0)+iRod);// has to seed every for loop otherwise the random numbers in different loops will be the same
//     rnu.seed((unsigned int)time(0)+iRod);// has to seed every for loop otherwise the random numbers in different loops will be the same
// 
//     double dtheta_avg = 0.0;
//     double dtheta2_avg = 0.0;
//     
//     blitz::Array<double,1> theta(nNodesperRod);
// 
// 
//     theta(0) = rnu.random();
//     theta(0) *= M_PI/2.0;
//     std::cout << "a = " << 0 
// 		<< " dtheta = " << 0.0
// 		<< " theta = " << theta(0)*180.0/M_PI 
// 		<< std::endl;
//     for ( int a=1; a<nNodesperRod; a++) {      
//       double dtheta = rng.random(); 
//       theta(a) = theta(a-1) + dtheta;
//       dtheta_avg += dtheta;
//       dtheta2_avg += sqr(dtheta);
//       std::cout << "a = " << a 
// 		<< " dtheta = " << dtheta*180.0/M_PI 
// 		<< " theta = " << theta(a)*180.0/M_PI 
// 		<< std::endl;
//     }
//     //theta(nNodesperRod-1) = 0.0;
//     
//     // make the rod straight
//     //theta = 0.0;
// 
// 
//     // temporary container for the nodes of a particular filament
//     SemiflexibleGel<2>::DefNodeContainer rod_nodes;
// 
//     // vector for nodal coordinates, initialized randomly for first node
//     BrownianNode<2>::Point X(0.0);
//     X = rnu.random()*Lbox, rnu.random()*Lbox;
//     for ( int a=0; a<nNodesperRod; a++) {
//       unsigned int id = a+nNodesperRod*iRod;
//       NodeBase::DofIndexMap idx(2);
//       idx[0]=2*id; idx[1]=2*id+1;
//       SemiflexibleGel<2>::DefNode * nd 
// 	= new SemiflexibleGel<2>::DefNode(id,idx,X,X);
//       
//       // add to this filament's node list
//       rod_nodes.push_back( nd );
//       // add to global node list
//       nodes.push_back( nd );
//       
//       BrownianNode<2>::Point dX(0.0);
//       dX = dL*cos(theta(a)), dL*sin(theta(a));
//       X = X + dX;
//     }
//     dtheta_avg /= nNodesperRod;
//     dtheta2_avg /= nNodesperRod;
//     
//     std::cout << " dtheta_avg = " << dtheta_avg 
// 	      << " dtheta2_avg = " << dtheta2_avg << std::endl;
//     
//     // add the new filament to the gel body
//     gel->addFilament( rod_nodes, kAngle, viscosity, kT, dt, kC, fitOrder);
// 
// 
//   }  // end of loop over filaments
// 
//   // Define periodic BC
//   LeesEdwards *box = new LeesEdwards (Lbox, Lbox, 0.25);
//  
//   gel->setBox( box ); 


  //////////////////////////////////////////////////////////////////////////
  //
  // Generate crosslinks, adding one wherever two filaments intersect.
  //
  //////////////////////////////////////////////////////////////////////////

//   const SemiflexibleGel<2>::FilamentContainer & filaments = gel->filaments();
//   SemiflexibleGel<2>::ConstFilamentIterator f1 = filaments.begin();
//   for(; f1!=filaments.end(); f1++){
//     SemiflexibleGel<2>::ConstFilamentIterator f2 = f1+1;
//     for(; f2!=filaments.end(); f2++){
//       SemiflexibleGel<2>::RodContainer rods1 = (*f1)->rods;
//       SemiflexibleGel<2>::RodContainer rods2 = (*f2)->rods;
//       SemiflexibleGel<2>::RodIterator r1 = rods1.begin();
//       SemiflexibleGel<2>::RodIterator r2 = rods2.begin();
//       for (; r1!=rods1.end(); r1++){
// 	r2 = rods2.begin();
// 	for (; r2!=rods2.end();r2++){
// 	  SemiflexibleGel<2>::DefNodeContainer rodnodes1 = (*r1)->getNodes();
// 	  SemiflexibleGel<2>::DefNodeContainer rodnodes2 = (*r2)->getNodes();
// 	  const BrownianNode<2>::Point & X1A = rodnodes1[0]->position();
// 	  const BrownianNode<2>::Point & X1B = rodnodes1[1]->position();
// 	  const BrownianNode<2>::Point & X2A = rodnodes2[0]->position();
// 	  const BrownianNode<2>::Point & X2B = rodnodes2[1]->position();
// 	  
// 	  Vector2D dx(0.0);
// 	  dx = 0.5*(X1A+X1B) - 0.5*(X2A+X2B);
// 	  box->mapDistance(dx);
// 	  if(norm2(dx) < 0.5*dL){
// 	    double xi1 = 0.5;
// 	    double xi2 = 0.5;
// 	    double k_crosslink = kBond;//0.0;
// 	    Crosslink<2> * crosslink 
// 	      = new Crosslink<2> ( rodnodes1[0], rodnodes1[1], 
// 				   rodnodes2[0], rodnodes2[1], 
// 				   xi1, xi2, 
// 				   k_crosslink, 
// 				   box );
// 	    gel->addCrosslink(crosslink);
//         std::cout << "crosslink created " << " X1A= " << X1A  
// 			<< " X1B= " << X1B << std:: endl
//                                           << " X2A= " << X2A  
//                                           << " X2B= " << X2B << std::endl;  
// 	  }
// 	}
//       }
//     }
//   }

  /////////////////////////////////////////////////////////////////////////////
  //
  // Create a model and a minimization solver and point them to the
  // SemiflexibleGel body.
  //
  /////////////////////////////////////////////////////////////////////////////

  Model::NodeContainer modelnodes;
  modelnodes.reserve( nodes.size() );
  for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
    modelnodes.push_back( *n );
  }
  Model model(modelnodes);
  model.pushBackBody( gel );


  //
  // initialize Quasi-Newton BFGS solver
  //
//   LeesEdwards * lebox = new LeesEdwards(syssize[0],syssize[1],shr);
//   gel->setBox(lebox);
//   if(kcl < 0.0) {
//     SemiflexibleGel<2>::ConstraintContainer ties = gel->constraints();
//     for(SemiflexibleGel<2>::ConstConstraintIterator ci = ties.begin(); ci != ties.end(); ci++) {
//       PeriodicTie<2>* pt = *ci;
//       pt->setBox(lebox);
//     }
//   }
//   else {
//     SemiflexibleGel<2>::CrosslinkContainer clks = gel->crosslinks();
//     for(SemiflexibleGel<2>::ConstCrosslinkIterator cl = clks.begin(); cl != clks.end(); cl++) {
//       (*cl)->setBox(lebox);
//     }
//   }
  int m=7;
  double factr=1.0e+4;
  double pgtol=1.0e-6;
  int iprint = 0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m ;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  int maxiter = 100;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxiter );
  
  GelOutput<2> output;
  output( gel, "relaxed-0" );

//   gel->print("relaxed-0");
  for(int step = 1; step<20; step++) {
    solver.solve( &model );
    char fname[100];
    sprintf(fname,"relaxed-%d",step);
//     gel->print(fname);
    output(gel, fname);
  }
  
    std::ofstream outStream(shrfilename,ios::app);
    outStream << shr << "\t" << gel->energy() << std::endl;
    outStream.close();
    delete gel;
    }
  std::cout << "All done." << std::endl;
  return 0;
  
}
