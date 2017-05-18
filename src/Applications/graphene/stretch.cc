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
  double extension_inp=0.0;
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

  double orientation=0;

  // afm radius
  double afmR = 175.0; // nm

  optind=2; // start at argv[2] looking for command-line options
  while( (option=getopt(argc,argv,"f:e:m:o:s:v:u:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 'e':
      extension_inp = std::atof(optarg);
      std::cout << "max extension: " << extension_inp << std::endl;
      break;
    case 'm':
      minStep_inp = std::atof(optarg);
      std::cout << "min step: " << minStep_inp << std::endl;
      break;
    case 'o':
      orientation = std::atof(optarg);
      std::cout << "orientation: " << orientation << std::endl;
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
  double scale = 1.0e3; // 1 micron = 1.0e3 nm
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;

  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;
    mesh->GetPoint(a, &(x[0]));
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  cout << "Number of nodes: " <<nodes.size() << endl;

  double C0 = 0.0;

  // find min and max x values
  double xmax=defNodes[0]->getPoint(0);
  double xmin=xmax;
  for(int i=1; i<defNodes.size(); i++) {
    double x=defNodes[i]->getPoint(0);
    xmax = std::max(x, xmax);
    xmin = std::min(x, xmin);
  }

  scale /=(xmax-xmin);
  // rescale size 
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = defNodes[i]->point();
    x *= scale;
    defNodes[i]->setPoint(x);
    defNodes[i]->setPosition(x);
  }

  xmin *= scale; 
  xmax *= scale;

  std::cout << "(xmin, xmax) = (" << xmin << "," << xmax << ")" << std::endl;

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
  Graphene stretching(orientation);

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

  // identify nodes on left and right edges
  std::vector<int> leftNodes;
  std::vector<int> rightNodes;
  for( int a=0; a<defNodes.size(); a++) {
    double L = xmax - xmin;
    double x = defNodes[a]->getPoint(0);

    if( x < xmin+0.005*L ) leftNodes.push_back(a);

    else if( x > xmax-0.005*L ) rightNodes.push_back(a);

  }

  // pin left nodes
  std::cout << "Left nodes: " <<std::endl;
  for( int a=0; a<leftNodes.size(); a++) {
    int A=leftNodes[a];    
    std::cout << A << std::endl;
    for(int i=0; i<3; i++) {
      nbd(3*A+i) = 2;
      hi(3*A+i) = lo(3*A+i) =  defNodes[A]->getPoint(i);
    }
  }

  // put right nodes on rollers
  std::cout << "Right nodes: " <<std::endl;
  for( int a=0; a<rightNodes.size(); a++) {
    int A=rightNodes[a];
    std::cout << A << std::endl;
    // x-dir
    nbd(3*A+0) = 2;
    hi(3*A+0) = lo(3*A+0) =  defNodes[A]->getPoint(0);
    // z-dir
    nbd(3*A+2) = 2;
    hi(3*A+2) = lo(3*A+2) =  defNodes[A]->getPoint(2);
  }

  solver.setBounds(nbd,lo,hi);

  std::cout << "Number of pin supports on left: " << leftNodes.size() 
	    << std::endl;

  std::cout << "Number of roller supports on left: " << rightNodes.size() 
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


  // add some viscosity for regularization
  //
  // f_v = \mu * dx
  //
  // set viscosity relative to Augmented Lagrangian contact stiffness
  double viscosity = viscosity_inp;
  std::cout << "Viscosity = " << viscosity << std::endl;
  ViscousRegularizer vr(bdc[0]->nodes(), viscosity);
  bdc[0]->addElement( &vr ); 

  double vrTol = 1.0e-6;

  double vrEnergy = vr.energy();
  double bdEnergy = bdc[0]->energy();

  string fzName = modelName + ".fz";
  // ofstream FvsZ(fzName.c_str());

  double F_prev = -1.0;
  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
  u_prev = 0.0;

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%

  int step=0;
  double dx = step_inp;
  for(double x = xmax; x<xmax+extension_inp; x+=dx, step++) {

    // move right nodes on rollers
    for( int a=0; a<rightNodes.size(); a++) {
      int A=rightNodes[a];
      hi(3*A+0) = lo(3*A+0) =  x;
    }
    solver.setBounds(nbd,lo,hi);
    
    // initial guess
    if( step > 0 ){ // 
      // add scaled version of previous displacement as an initial
      // guess, but only in loading direction
      for(int i=0; i<model.dof(); i++ ) {
	solver.field(i) = x_prev(i) + 0.99*dx*u_prev(i);
      }      
      model.putField(solver);
    }

    // step forward in "time", relaxing viscous energy & forces 
    vr.step();
    
    std::cout << std::endl
	      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	      << std::endl
	      << " x = " << x 
	      << " extension = " << x-xmax 
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


      model.computeAndAssemble(solver, false, true, false);
      
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


      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < std::abs(vrTol*bdEnergy) && 
	 solver.projectedGradientNorm()<=pgtol ) {
	// viscous energy is small enough; exit
	break;
      }
    }

    if(verbose) {
      std::cout << "Contact and Viscous energy converged." << std::endl
		<< std::endl
		<< "   extension = " << x-xmax << std::endl
		<< std::endl;
      }


    // Compute displacement
    for(int i=0; i<model.dof(); i++ ) {
      u_prev(i) = ( solver.field(i) - x_prev(i) )/std::abs(dx);
    }

    // Save current (successful) state as previous
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);

    for(int b=0; b<bdc.size(); b++) {
      char name[100]; 
      sprintf(name,"%s-step%04d",modelName.c_str(),step);
      bdc[b]->printParaview(name);
    }
    // FvsZ << std::setw( 24 ) << std::setprecision(16) 
    // 	 << Z
    // 	 << std::setw( 24 ) << std::setprecision(16) 
    // 	 << afm->FZ()
    // 	 << std::setw( 10 ) 
    // 	 << step
    // 	 << std::endl;

    // // check if we are done
    // if( unload && Z+dZ > Zbegin ) 
    //   break;
    // else if( !unload && Z+dZ < Zend ) 
    //   break;
  }
  // FvsZ.close();

  std::cout << "Indentation complete." << std::endl;

  return 0;
}

