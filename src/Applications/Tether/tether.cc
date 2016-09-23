#include <string>
#include <iostream>
#include <vector>
#include "Node.h"
#include "EvansElastic.h"
#include "FVK.h"
#include <tvmet/Vector.h>
#include <fstream>
#include "LoopShellBody.h"
#include "Model.h"
#include "Solver.h"
#include "ConjugateGradientWSK.h"
#include "ViscousRelaxation.h"
#include "SimulatedAnnealing.h"
#include "Rigid.h"

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);


int main(int argc, char* argv[])
{
  bool verbose=true;

#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif


  ifstream ifs;
  string ofn;
 
  ioSetting(argc, argv, ifs, ofn);
	
  //
  // create vector of nodes
  int dof=0;
  std::vector< NodeBase* > nodes;
  char key;
  ifs>>key;
  DeformationNode<3>::Point xc;
  xc = 0.0, 0.0, 0.0;  
  double R=0.0;
  for ( int i = 0; key=='v'; i++, ifs>>key){
    int id=i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    xc += x;
    R = std::max(R,x(0));
    R = std::max(R,x(1));

    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    nodes.push_back(new DeformationNode<3>(id,idx,x));
  }
  xc /= nodes.size();

//   //
//   // make vesicle spherical
//   //
//   for(std::vector<NodeBase*>::iterator n=nodes.begin(); n!=nodes.end(); n++){
//     double r1 = (*n)->getPoint(0) - xc(0);
//     double r2 = (*n)->getPoint(1) - xc(1);
//     double r3 = (*n)->getPoint(2) - xc(2);
//     double x3 = xc(2);
//     if(r3 < 0.0) x3 -= sqrt(std::max(0.0,R*R - r1*r1 - r2*r2));
//     else x3 += sqrt(R*R - r1*r1 - r2*r2);
//     (*n)->setPoint(2,x3);
//     static_cast<DeformationNode<3>*>(*n)->resetPosition();
//   }

  double Min=std::numeric_limits<double>::max();
  double Max=-std::numeric_limits<double>::max();
  if(verbose) cout << nodes.size() << endl;

  for(int i=0; i<nodes.size(); i++) {
    Min = std::min(Min,nodes[i]->getPoint(2));
    Max = std::max(Max,nodes[i]->getPoint(2));
  }
  // Set Zmin and Zmax just inside of the actual ends
  double Zmin = Min + 0.10*(Max-Min);
  double Zmax = Max - 0.10*(Max-Min);

//   // set Zmin and Zmax to heights of top and bottom 1-rings
//   double Zmin=std::numeric_limits<double>::max();
//   double Zmax=-std::numeric_limits<double>::max();
//   for(int i=0; i<nodes.size(); i++) {
//     double Z = nodes[i]->getPoint(2);
//     if( Z < Max && Z > Min ) {
//       Zmin = std::min(Zmin,Z);
//       Zmax = std::max(Zmax,Z);
//     }
//   }

  //
  // create connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  for (int i = 0; key=='f'; i++){
    int tmp=0;
    if((string)(argv[3])=="neg") {
      ifs >> tmp; c[1]=tmp-1;
      ifs >> tmp; c[0]=tmp-1;
      ifs >> tmp; c[2]=tmp-1;
    } else {
      ifs >> tmp; c[0]=tmp-1;
      ifs >> tmp; c[1]=tmp-1;
      ifs >> tmp; c[2]=tmp-1;  
    }
    connectivities.push_back(c);
    if(!(ifs>>key)) break;
  }
  if(verbose) cout << connectivities.size() << endl;
  ifs.close();
	
  double KC = 1.0e0;
  double KG = 0.0;
  double C0 = 0.0;
  double mu = 1.0e0;
  double KS = 1.0e0;

  double viscosity=1.0e-2;

  ifstream inp("parameters.inp");
  inp >> KC >> KG >> mu >> KS >> viscosity;
  if(verbose) {
    std::cout << "Input KC: " << KC << std::endl
	      << "Input KG: " << KG << std::endl
	      << "Input mu: " << mu << std::endl
	      << "Input KS: " << KS << std::endl
	      << "Input viscosity: " << viscosity << std::endl;
  }
  typedef EvansElastic MaterialType;
  MaterialType bilayer( KC, KG, C0, mu, KS );
//   typedef FVK MaterialType;
//   MaterialType bilayer( KC, KG, C0, mu, 0.3 );
  if(verbose) 
    std::cout << "SCElastic Material has been created." << std::endl;
	
  //
  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  double penaltyVolume=1.0e7;
  double penaltyArea=1.0e7;
  int quadOrder = 2;
  
  typedef LoopShellBody<MaterialType> LSB;
  LSB bd(bilayer, connectivities, nodes, quadOrder, 
	 nBoundaries, pressure, tension,
	 penaltyVolume, penaltyArea, viscosity,
	 LoopShell<MaterialType>::penalty, 
	 LoopShell<MaterialType>::penalty 
	 );
  if(verbose)
    cout << "Created a body." << endl
	 << " volume = " << bd.volume() << endl;

  bd.setOutput(Body::paraview);

  bd.compute(true,true,false);

  double k2 = 1.0e2;
  inp >> k2;
  std::vector< PenaltyBC* > topBCs, botBCs;
  for(int i=0; i<nodes.size(); i++) {
    DeformationNode<3> * apex = dynamic_cast< DeformationNode<3>* >(nodes[i]);
    if( apex != 0 ) {
      double X = apex->getPoint(0);
      double Y = apex->getPoint(1);
      double Z = apex->getPoint(2);
      if( ( Z >= Zmax || Z <= Zmin ) /*&& X == 0.0 && Y == 0.0*/ ) {
	if(verbose)
	  std::cout << "Fixing node " << apex->id()
		    << " at  " << X << ", " << Y << ", " << Z << std::endl;
	PenaltyBC::Vector3D x0(0.0);
	x0 = X,Y,Z;
	PenaltyBC::VectorBC bc; bc = true, true, true;
	PenaltyBC * apexBC = 
	  new PenaltyBC( apex, x0, bc, k2 );
	bd.pushBack( apexBC );
	if( Z>=Zmax ) topBCs.push_back( apexBC );
	else 	      botBCs.push_back( apexBC );
      }
    }
  }
  if( topBCs.size() == 0 || botBCs.size() == 0 ) {
    std::cout << "topBCs.size() = " << topBCs.size() << std::endl;
    std::cout << "botBCs.size() = " << botBCs.size() << std::endl;
    return 0;
  }
  
  //
  // create Body Container
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  //
  // create Model
  Model model(bdc);

#ifdef WITH_MPI
  if(procId==0) 
#endif
    {
      string fname = argv[1];
      fname += ".shape";
      model.print(fname);
    }
//   model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
// #ifdef WITH_MPI
//   MPI_Finalize();
// #endif
//   return 0;
	
  ConjugateGradientWSK CGsolver;//(true);
  int maxIter = 100*model.dof();
  int restartStride = 2*model.dof();
  int printStride = 100;
  double tol = 1.0e-8;//8;
  double absTol = 1.0e-8;//8;
  double tolLS = 1.0e-6;
  int maxIterLS = 20;

  ifstream solverinp("solver.inp");
  solverinp >> maxIter >> printStride >> tol >> absTol;
  if(verbose)
    std::cout << "maxIter: " << maxIter << std::endl
	      << "printStride: " << printStride << std::endl
	      << "tol: " << tol << std::endl
	      << "absTol: " << absTol << std::endl;
   
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 maxIter, restartStride, printStride, 
			 tol, absTol, tolLS, maxIterLS);
  CGsolver.setWolfeParameters(1.0e-6, 1.0e-3, 1.0e-4);

  model.print("Before");
  
  double dt=1.0e-8;
  solverinp >> dt;
  if(verbose) std::cout << "dt: " << dt << std::endl;
  ViscousRelaxation VRsolver(dt,tol,absTol,maxIter,printStride);

  Solver * solver;
  if( (string)(argv[2])=="ev" ) {
    model.checkRank(model.dof(),true);
    return 0;
  } else if( (string)(argv[2])=="vr" ) {
    solver = &VRsolver;
  } else {
    solver = &CGsolver;
  }

  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = bd.constraintVolume(); 
  double A = bd.constraintArea();  
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);

  double tolV = 1.0e-6;
  inp >> tolV;
  if(verbose) std:: cout << "tolV: " << tolV << std::endl;
  int maxIterV = 100;
  inp >> maxIterV;
  if(verbose) std:: cout << "maxIterV: " << maxIterV << std::endl;

  double Vfactor=1.0;
  double Afactor=1.0;
  inp >> Vfactor >> Afactor;
  if(verbose) 
    std::cout << "Vfactor: " << Vfactor << std::endl
	      << "Afactor: " << Afactor << std::endl;

  bd.setConstraintVolume(Vfactor*V);
  bd.setConstraintArea(Afactor*A);

  bd.compute(true,false,false);
  if(verbose) 
    std::cout << std::setprecision(16) 
	      << "Before:"<< std::endl
	      << "Volume = "<< bd.volume() << std::endl
	      << "Ref. Volume = "<< bd.constraintVolume() << std::endl
	      << "Area = "<< bd.area() << std::endl
	      << "Ref. Area = "<< bd.constraintArea() << std::endl
	      << "Reduced Volume = " << vred << std::endl
	      << "Ref. Reduced Volume = " << Vred << std::endl
	      << "Energy = " << bd.energy() << std::endl;
  
  for(int iter=0; iter<maxIterV; iter++) {
    if(verbose) std::cout << "iter = "<< iter << std::endl;
    bd.resetReference();
//     CGsolver.solve( &model );
    solver->solve( &model );
    double Ftop = 0.0;
    for(std::vector< PenaltyBC* >::const_iterator bc=topBCs.begin();
	bc != topBCs.end(); bc++ ) 
      Ftop += (*bc)->force(2);
    double Fbot = 0.0;
    for(std::vector< PenaltyBC* >::const_iterator bc=botBCs.begin();
	bc != botBCs.end(); bc++ ) 
      Fbot += (*bc)->force(2);
//     std::cout << "Ftop: " << Ftop << std::endl
//	<< "Fbot: " << Fbot << std::endl;

    double constraintEnergy = bd.constraintEnergy();
    double viscousEnergy = bd.viscousEnergy();
    double energy = bd.energy();

#ifdef WITH_MPI
    double tmp = Ftop;
    MPI_Allreduce(&tmp, &Ftop, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
    tmp = Fbot;
    MPI_Allreduce(&tmp, &Fbot, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
    tmp = constraintEnergy;
    MPI_Allreduce(&tmp, &constraintEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
    tmp = viscousEnergy;
    MPI_Allreduce(&tmp, &viscousEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
    tmp = energy;
    MPI_Allreduce(&tmp, &energy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
#endif

    if(verbose) 
      std::cout << "c: " << constraintEnergy << std::endl
		<< "v: " << viscousEnergy << std::endl
		<< "e: " << energy << std::endl
		<< "Top force    = " << Ftop << std::endl
		<< "Bottom force = " << Fbot << std::endl
		<< "Pressure = " << bd.pressure() << std::endl
		<< "Tension = " << bd.tension() << std::endl;

    if( constraintEnergy + viscousEnergy < tolV*energy ) 
      break;
  }
  
  CGsolver.setParameters(voom::ConjugateGradientWSK::Secant,
			 voom::ConjugateGradientWSK::PR,
			 model.dof(), restartStride, printStride, 
			 tol/100, absTol, tolLS, maxIterLS);
  solver->solve( &model );

  double Ftop = 0.0;
  for(std::vector< PenaltyBC* >::const_iterator bc=topBCs.begin();
      bc != topBCs.end(); bc++ ) 
    Ftop += (*bc)->force(2);
  double Fbot = 0.0;
  for(std::vector< PenaltyBC* >::const_iterator bc=botBCs.begin();
      bc != botBCs.end(); bc++ ) 
    Fbot += (*bc)->force(2);

#ifdef WITH_MPI
  double tmp = Ftop;
  MPI_Allreduce(&tmp, &Ftop, 1, MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);
  tmp = Fbot;
  MPI_Allreduce(&tmp, &Fbot, 1, MPI_DOUBLE, 
		MPI_SUM, MPI_COMM_WORLD);
  if(procId != 0) {
    MPI_Finalize();
    return 0;
  }
#endif

  std::string fname = argv[1];
  fname += ".relaxed";
  model.print(fname);
  bd.printObj(fname);
  v = bd.volume();
  a = bd.area(); 
  vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  V = bd.constraintVolume();
  A = bd.constraintArea(); 
  Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
  if(verbose) 
    std::cout   << "Volume = "<< bd.volume() << std::endl
		<< "Ref. Volume = "<< bd.constraintVolume() << std::endl
		<< "Area = "<< bd.area() << std::endl
		<< "Ref. Area = "<< bd.constraintArea() << std::endl
		<< "Reduced Volume = " << vred << std::endl
		<< "Ref. Reduced Volume = " << Vred << std::endl
		<< "Energy = " << solver->function() << std::endl
		<< "Top force    = " << Ftop << std::endl
		<< "Bottom force = " << Fbot << std::endl
		<< "Pressure = " << bd.pressure() << std::endl
		<< "Tension = " << bd.tension() << std::endl;
  fname = argv[1];
  fname += ".info";
  ofstream info(fname.c_str());
  info << "Volume = "<< bd.volume() << std::endl
       << "Ref. Volume = "<< bd.constraintVolume() << std::endl
       << "Area = "<< bd.area() << std::endl
       << "Ref. Area = "<< bd.constraintArea() << std::endl
       << "Reduced Volume = " << vred << std::endl
       << "Ref. Reduced Volume = " << Vred << std::endl
       << "Energy = " << solver->function() << std::endl
       << "Top force    = " << Ftop << std::endl
       << "Bottom force = " << Fbot << std::endl
       << "Pressure = " << bd.pressure() << std::endl
       << "Tension = " << bd.tension() << std::endl;

#ifndef WITH_MPI
  fname = argv[1];
  fname += ".FvsZ";
  ofstream FvsZ(fname.c_str());

  double Zlen = Zmax-Zmin;
  double Zstep = Zlen/100.0;
  for(double slice=Zmin; slice<=Zmax; slice+=Zstep) {
    double beadForce=0.0;
    for(LSB::ConstNodeIterator n=bd.nodes().begin(); n!=bd.nodes().end(); n++) {
      LSB::FeNode_t * sn = dynamic_cast<LSB::FeNode_t*>(*n);
      if( sn ) {
	if( sn->getPoint(2) >= slice ) beadForce += sn->getForce(2);
      }
    }

    FvsZ << std::setw(30) << slice 
	 << std::setw(30) << beadForce << std::endl;
  }
#endif

#ifdef WITH_MPI
  MPI_Finalize();
#endif
  return 0;
}


void ioSetting(int argc, char* argv[], ifstream& ifs, string& ofn)
{

  // if no arguement or too many
  if (argc < 2 || argc > 4)
    {
      cout << "Usage: ProgramName InputFileName [cg|vr] [pos|neg]" << endl;
      exit(0);
    }

  string inFullName = argv[1];
  string pathName;

  // create input stream
  ifs.open( inFullName.c_str(), ios::in);
  if (!ifs)
    {
      cout << "can not open input file: " << inFullName << endl;
      exit(0);
    }
	
}
