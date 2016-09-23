#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "C0MembraneBody.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Contact.h"
#include "RigidHemisphereAL.h"
#include "ViscousRegularizer.h"
#include "Graphene.h"

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;



int main(int argc, char* argv[])
{
  if( argc < 2 ) {
      cout << "Usage: indent modelName [-n nsteps]." << endl;
      return(0);
  }

#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  bool verbose=true;
#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif

  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

  string modelName = argv[1];

  double gamma_inp=0.0;
  double indent_inp=0.0;
  double viscosity_inp=0.0;
  double friction_inp=0.0;
  double minStep_inp=0.01;
  double step_inp=0.0;
  double prestrain=0.0;
  string prestressFlag = "yes";
  bool CST=false;
  bool badCommandLine=false;
  bool unload=false;
  int option;

  // afm radius
  double afmR = 175.0; // nm

  optind=2; // start at argv[2] looking for command-line options
  while( (option=getopt(argc,argv,"f:i:m:p:s:v:u:R:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 'i':
      indent_inp = std::atof(optarg);
      std::cout << "max indentation: " << indent_inp << std::endl;
      break;
    case 'm':
      minStep_inp = std::atof(optarg);
      std::cout << "min step: " << minStep_inp << std::endl;
      break;
    case 's':
      step_inp = std::atof(optarg);
      std::cout << "step size: " << step_inp << std::endl;
      break;
    case 'v':
      viscosity_inp = std::atof(optarg);
      std::cout << "viscosity: " << viscosity_inp << std::endl;
      break;
    case 'f':
      friction_inp = std::atof(optarg);
      std::cout << "friction: " << friction_inp << std::endl;
      break;
    case 'u':
      unload = false;
      std::cout << "Unloading curve will be calculated." << std::endl;
      break;
    case 'p':
      prestrain = std::atof(optarg);
      std::cout << "prestrain: " << prestrain << std::endl;
      break;
    case 'R':
      afmR = std::atof(optarg);
      std::cout << "afm radius: " << afmR << std::endl;
      break;
    case '?':
      std::cout << "Option " << char(optopt) << " is not supported." 
		<< std::endl;
      return 1;
    default :
      badCommandLine=true;
      break;
    }
  }

  std::cout << "optind = " << optind << std::endl
	    << "option = "<< option << std::endl;

  string inputFileName = modelName + ".vtk";
  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  // send through normals filter to ensure that triangle orientations
  // are consistent
  // vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  // normals->SetInput( reader->GetOutput() );
  // normals->ConsistencyOn();
  // normals->SplittingOff();
  // normals->AutoOrientNormalsOn();
  vtkDataSet * mesh = reader->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  // create vector of nodes
  double Rdisc = 5.0e2; // 1 micron = 1.0e3 nm
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;
  double Rmax = 0;

  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;
    mesh->GetPoint(a, &(x[0]));
    double r=tvmet::norm2(x);
    Ravg += r;
    Rmax = std::max(r,Rmax);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  cout << "Number of nodes: " <<nodes.size() << endl
       << "Ravg = " << Ravg << endl
       << "Rmax= " << Rmax << endl;


  double C0 = 0.0;

  // rescale size 
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = defNodes[i]->point();
    x *= Rdisc/Rmax;
    defNodes[i]->setPoint(x);
    x -= prestrain*x;
    defNodes[i]->setPosition(x);
  }

  double KC = 0.0;

  //
  // Cadelano, et al., PRL 102: 235502 (2009)
  //
  // double Y = 312.0e3; // 312 N/m = 312 nN/nm = 312e3 pN/nm
  // double nu = 0.31;

  // double c111=-2724.7e3;
  // double c222=-2523.2e3;
  // double c112=-591.1e3;

  //
  // Wei, et al., PRB 80: 205407 (2009)
  //
  // double Y = 348.0e3; 
  // double nu = 0.169;

  // double c111=-2817.0e3;
  // double c222=-2693.3e3;
  // double c112=-337.1e3;

  
  // create Model
  Model::BodyContainer bdc;


  // read in triangle connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=mesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  if(verbose) cout << "Number of triangles: " <<ntri << endl;

  int nodes_per_element=mesh->GetCell(0)->GetNumberOfPoints();
  assert( nodes_per_element == 3 || nodes_per_element == 6 );

  ShapeTri3 * shape3 = new ShapeTri3( Shape<2>::CoordinateArray(0.0) );
  ShapeTri6 * shape6 = new ShapeTri6( Shape<2>::CoordinateArray(0.0) );

  Shape<2> * shape=shape3;
  
  if (nodes_per_element == 6 ) {
    shape = shape6;
  }

  unsigned int quadOrder = 1;
  if( nodes_per_element == 6 ) quadOrder = 2;

  cout << "Creating TriangleQuadrature." << endl;

  // stretching body
  std::vector< std::vector<int> > s_connectivities;
  for (int e = 0; e<ntri; e++){
    std::vector<int> c(nodes_per_element);
    for(int a=0; a<nodes_per_element; a++) 
      c[a]= mesh->GetCell(e)->GetPointId(a);
    s_connectivities.push_back(c);
  }

  // MaterialType stretching( 0.0, 0.0, C0, Y, nu ); 
  // MaterialType stretching( Y, nu, c111, c222, c112 );	
  Graphene stretching;

  if( nodes_per_element == 3 ) {
    typedef C0MembraneBody<TriangleQuadrature,Graphene,ShapeTri3> MB;
    MB * bdm = new MB(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm->setOutput(paraview);
  
    Vector3D N;
    N = 1.0, 0.0, 0.0;
    bdm->addStressDirection(N);
    N = 0.0, 1.0, 0.0;
    bdm->addStressDirection(N);
    
    bdc.push_back(bdm);

  } else {
    typedef C0MembraneBody<TriangleQuadrature,Graphene,ShapeTri6> MB;
    MB * bdm = new MB(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm->setOutput(paraview);
  
    Vector3D N;
    N = 1.0, 0.0, 0.0;
    bdm->addStressDirection(N);
    N = 0.0, 1.0, 0.0;
    bdm->addStressDirection(N);
    
    bdc.push_back(bdm);
  }


  // bdc[0]->checkConsistency();
  // return 0;

  Model model(bdc,nodes);


  int m=5;
  int maxIter=500;//1000;//model.dof();
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = 0;
  int maxIter_inp=0;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m >> maxIter_inp;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input maxIter: " << maxIter_inp << std::endl;
  maxIter = std::max(maxIter,maxIter_inp);

  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );//(true);

  //
  // Here.  Next steps:
  // 1. fix nodes on the boundary using solver bounds
  // 2. clean up rest of the file so it compiles
  // 3. look up tip size and indentation amplitude from emails with Haider
  //

  // set up bounds for solver
  bool boundaryConditionFlag = true;
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;
  int nFixedNodes=0;
  for( int a=0; a<defNodes.size(); a++) {
    if( tvmet::norm2( defNodes[a]->point() ) < Rdisc*0.99 ) continue;
 
    nFixedNodes++;

    for(int i=0; i<3; i++) {
      nbd(3*a+i) = 2;
      hi(3*a+i) = lo(3*a+i) =  defNodes[a]->getPoint(i);
    }
  }

  solver.setBounds(nbd,lo,hi);

  std::cout << "Number of fixed nodes on boundary: " << nFixedNodes 
	    << std::endl;

  for(int n=0; n<nodes.size(); n++) {
    for(int i=0; i<nodes[n]->dof(); i++) nodes[n]->setForce(i,0.0);
  }
    
  for(int b=0; b<bdc.size(); b++) {
    std::cout << "bdc[" << b << "]->compute()" << std::endl;
    bdc[b]->compute(true,true,false);    
  }

  std::cout << "Initial Shape." << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  string fname = modelName;
  fname += ".initial";
  model.print(fname);

  // relax initial shape
  solver.solve(&model);

  std::cout << "Shape relaxed." << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  fname = modelName;
  fname += ".relaxed";
  model.print(fname);


  // create indentor and plate

  tvmet::Vector<double,3> xc; xc = 0.0, 0.0, afmR; 

  double friction = friction_inp;

//   RigidHemisphereContact * afm 
//     = new RigidHemisphereContact( defNodes, afmR, xc, friction, pgtol );
//   bd->pushBackConstraint( afm );
//   RigidHemisphereElement * afm 
//     = new RigidHemisphereElement( defNodes, afmR, xc, friction, pgtol );
  double k_AL = 1.0e4;
  RigidHemisphereAL * afm 
    = new RigidHemisphereAL( defNodes, k_AL, afmR, xc, friction );
  afm->updateContact();
  
  bdc[0]->addElement( afm );
  std::cout << "Added afm to body." << std::endl;


  double Zbegin = 0.0;
  double Zend = Zbegin-std::abs(indent_inp);
  double dZ = (Zend-Zbegin)/100;
  if(step_inp > 0.0) {
    dZ = -std::abs(step_inp);
  }

  std::cout << "Zbegin = " << Zbegin << std::endl
	    << "Zend = " << Zend << std::endl
	    << "dZ = " << dZ << std::endl;

  // add some viscosity for regularization
  //
  // f_v = \mu * dx
  //
  // set viscosity relative to Augmented Lagrangian contact stiffness
  double viscosity = viscosity_inp*k_AL;
  std::cout << "Viscosity = " << viscosity << std::endl;
  ViscousRegularizer vr(bdc[0]->nodes(), viscosity);
  bdc[0]->addElement( &vr ); 

  double vrTol = 1.0e-6;

  double vrEnergy = vr.energy();
  double bdEnergy = bdc[0]->energy();

  string fzName = modelName + ".fz";
  ofstream FvsZ(fzName.c_str());

  double F_prev = -1.0;
  double Z_drop = Zbegin;
  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
  u_prev = 0.0;

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%

  int step=0;
  for(double Z = Zbegin; /*Z<Zend+0.5*dZ*//*Z>=Zbegin*/; Z+=dZ, step++) {
    
    // initial guess
    if( step > 0 ){ // 
      // add scaled version of previous displacement as an initial
      // guess, but only in loading direction
      for(int i=0; i<model.dof(); i++ ) {
	solver.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
      }      
      model.putField(solver);
    }

    // step forward in "time", relaxing viscous energy & forces 
    vr.step();
    

    // move afm down by dZ
    afm->setZ(Z+afmR);


    std::cout << std::endl
	      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	      << std::endl
	      << " Z = " << Z 
	      << " zeta = " << Z-Zbegin 
	      << std::endl;

    //viscosity = std::max(viscosity,viscosity_inp);
    //vr.setViscosity(viscosity_inp);
    int viterMax = 20;
    for(int viter=0; viter<viterMax; viter++) {

      blitz::Array<double,1> vSave(model.dof());
      model.getField(solver);
      for(int i=0; i<model.dof(); i++ ) vSave(i) = solver.field(i);
      
      if(verbose) std::cout << std::endl 
			    << "VISCOUS ITERATION: " << viter 
			    << "\t viscosity = " << vr.viscosity()
			    << std::endl
			    << std::endl;


//       if(viter == 0) {
// 	glass->updateContact();
// 	afm->updateContact();	
// 	model.computeAndAssemble(solver, false, true, false);	
//       }

      // update contact
      afm->updateContact();


      model.computeAndAssemble(solver, false, true, false);
//      model.print("contact");

      
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bdc[0]->energy();

      if(verbose) {
	std::cout << "ENERGY:" << std::endl
		  << "viscous energy = " << vrEnergy << std::endl
		  << "   body energy = " << bdEnergy << std::endl
		  << "  total energy = " << solver.function() << std::endl
		  << std::endl;
      }

      // update contact
      if(verbose) {
	std::cout << "CONTACT:" << std::endl;
	std::cout << "       top active  = " << afm->active() << std::endl
		  << "        top force  = " << afm->FZ() << std::endl 
 		  << "  top penetration  = " << afm->penetration() << std::endl;
      }
	
      // update viscosity
//       if( vr.velocity() > 2*targetVelocity && vr.viscosity() < maxViscosity ) {

// 	// Displacements very large; re-do last step with max viscosity
// 	std::cout << std::endl
// 		  << "Velocity too large.  Re-do previous step with max viscosity." 
// 		  << std::endl << std::endl;
// 	for(int i=0; i<model.dof(); i++ ) solver.field(i) = vSave(i);
// 	model.putField(solver);
// 	vr.setViscosity(maxViscosity);	
//       } else {
// 	// adjust viscosity to get velocity equal to target
// 	double viscosity = vr.viscosity()*vr.velocity()/targetVelocity;
// 	viscosity = std::max(viscosity,minViscosity);
// 	vr.setViscosity(viscosity);
//       }

      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < std::abs(vrTol*bdEnergy) && solver.projectedGradientNorm()<=pgtol &&
	 afm->penetration() < 1.0e-2*std::abs(dZ) ) {
	// viscous energy is small enough; exit
	break;
      }
    }

    if(verbose) {
      std::cout << "Contact and Viscous energy converged." << std::endl
		<< std::endl
		<< "   indentation = " << Z << std::endl
		<< std::endl;
      }


    // Compute displacement
    for(int i=0; i<model.dof(); i++ ) {
      u_prev(i) = ( solver.field(i) - x_prev(i) )/std::abs(dZ);
    }

    // Save current (successful) state as previous
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);

    F_prev = afm->FZ();


    for(int b=0; b<bdc.size(); b++) {
      char name[100]; 
      sprintf(name,"%s-step%04d",modelName.c_str(),step);
      bdc[b]->printParaview(name);
    }
    FvsZ << std::setw( 24 ) << std::setprecision(16) 
	 << Z
	 << std::setw( 24 ) << std::setprecision(16) 
	 << afm->FZ()
	 << std::setw( 10 ) 
	 << step
	 << std::endl;

    // check if we are done
    if( unload && Z+dZ > Zbegin ) 
      break;
    else if( !unload && Z+dZ < Zend ) 
      break;
  }
  FvsZ.close();

  std::cout << "Indentation complete." << std::endl;

  return 0;
}

