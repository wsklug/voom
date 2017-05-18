#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShellPreconditioner.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "ViscousRegularizer.h"
#include "Spring.h"

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkWarpVector.h>
#include <vtkLoopSubdivisionFilter.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  if( argc < 3 ) {
      cout << "Usage: relax modelName FvK [faces]." << endl
	   << endl
	   << "\t [faces] is an optional list of integers between 0 and 11"
	   << endl
	   << "\t identifying faces of the dodecahedron to be removed."
	   << endl << endl;
	   
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

  // read in vtk file
  string modelName = argv[1];

  string inputFileName = modelName + ".vtk";

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  // send through normals filter to ensure that triangle orientations
  // are consistent
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  normals->SetInput( reader->GetOutput() );
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  vtkPolyData * mesh = normals->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  // create vector of nodes, get positions from vtk input
  double R = 1.0;
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;

  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;
    mesh->GetPoint(a, &(x[0]));
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  if(verbose) cout << "Number of nodes: " <<nodes.size() << endl
		   << "Ravg = " << Ravg << endl;


  // read in connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=mesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  if(verbose) cout << "Number of triangles: " <<ntri << endl;

  ////////////////////////////////////////////////////////////////////
  // Skip elements on some of the dodecahedral faces.
  // 
  // 1. Check to make sure that the mesh appears to be of a
  // dodecahedron.  Number of elements should be an integer multiple
  // of 12, i.e., nElements=n*12.
  //
  int nFaces=12;
  float triPerFace=0.0;
  float remainder=0.0;
  remainder = modf( ntri/float(nFaces), &triPerFace );
  if( remainder !=0.0 ) {
    std::cout << "ntri/12 = " << ntri/float(nFaces) << " is not an integer."
	      << std::endl;
    return 0;
  }
  // 
  // 2. Deactivate the first n=nElements/12 elements.
  //
  vector<bool> removed(nFaces,false);
  // to remove face N add below
  //      removed[N] = true;
  for(int arg=3; arg<argc; arg++) {
    int faceId = atoi(argv[arg]);
    removed[faceId]=true;
  }

  for(int face=0; face<12; face++) {
    if( removed[face] ) continue;
    // loop through triangles on each pentamer
    for( int e=0; e<triPerFace; e++ ) {
      int tri = triPerFace*face + e;
      assert(mesh->GetCell(tri)->GetNumberOfPoints() == 3);
      for(int a=0; a<3; a++) c[a] = mesh->GetCell(tri)->GetPointId(a);
      connectivities.push_back(c);
    }
  }

  // rescale reference such that average radius is R
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = defNodes[i]->position();
    x *= R/Ravg;
    defNodes[i]->setPosition(x);
//     x = defNodes[i]->point();
//     x *= R/tvmet::norm2(x);
    defNodes[i]->setPoint(x);
  }


  ////////////////////////////////////////////////////////////////////
  // Solution section
  ////////////////////////////////////////////////////////////////////


  // loop over a range of FvK numbers, evenly spaced on the log axis
  double gamma = atof(argv[2]);

  // initialize material response
  double scalingFactor = 1.0e6;

  double Y = sqrt(gamma);
  double KC = 1.0/Y;
  Y  *= scalingFactor;
  KC *= scalingFactor;
  double nu = 1.0/3.0;
  
  double KG =-KC ;
  double C0 = 0.0;
  
  typedef FVK MaterialType;
  
  bool springsON=false;
  
  double Yb=0.0e-4*Y;//Y;
  if(springsON) Yb=0.0;
  MaterialType bending( KC, KG, C0, Yb, nu );
  
  // create Body for bending
  unsigned int quadOrder = 1;

  typedef LoopShellBody<MaterialType> LSB;
  LSB bd(bending, connectivities, nodes, quadOrder);

  bd.setOutput(paraview);

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);
  
  if(springsON) {
    // create springs
    double d0=0.0;
    double dx;
    dx = norm2(   defNodes[connectivities[0](0)]->point() 
		- defNodes[connectivities[0](1)]->point()  );
    d0 = std::max(d0, dx);
    dx = norm2(   defNodes[connectivities[0](1)]->point() 
	        - defNodes[connectivities[0](2)]->point()  );
    d0 = std::max(d0, dx);
    dx = norm2(   defNodes[connectivities[0](2)]->point() 
		- defNodes[connectivities[0](0)]->point()  );
    d0 = std::max(d0, dx);

    std::cout << "d0 = " << d0 << std::endl;

    double k=0.25*sqrt(3.0)*Y;
    for(int f=0; f<connectivities.size(); f++) {
      DeformationNode<3> * nA = defNodes[connectivities[f](0)];
      DeformationNode<3> * nB = defNodes[connectivities[f](1)];
      DeformationNode<3> * nC = defNodes[connectivities[f](2)];
      bd.pushBack( new Spring<3>(nA,nB,k) );
      bd.pushBack( new Spring<3>(nB,nC,k) );
      bd.pushBack( new Spring<3>(nC,nA,k) );
    }
  } else {
    // create Body for stretching
    std::vector< std::vector<int> > s_connectivities;
    for(int i=0; i<connectivities.size(); i++) {
      std::vector<int> c(3);
      for(int j=0; j<3; j++) c[j]=connectivities[i](j);
      s_connectivities.push_back(c);
    }
    MaterialType stretching( 0.0, 0.0, C0, Y, nu );
    quadOrder = 1;
    typedef C0MembraneBody<TriangleQuadrature,MaterialType,ShapeTri3> MB;
    MB * bdm = new MB(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm->setOutput(paraview);
    bdc.push_back(bdm);
  } 

  ShellPreconditioner * pre = 
    new ShellPreconditioner(connectivities, defNodes, Y/10.0);

//   bdc.push_back(pre);

  Model model(bdc,nodes);

  // initialize Quasi-Newton BFGS solver
  int m=5;
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = -1;
  double pentol=1.0e-4;
  int maxIter = 2*model.dof();
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl;
  
#if 1
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );

    // set up bounds for solver
    bool boundaryConditionFlag = true;
    blitz::Array<int,1> nbd(3*nodes.size());
    blitz::Array<double,1> lo(3*nodes.size());
    blitz::Array<double,1> hi(3*nodes.size());
    nbd = 0;
    lo = -1.5;
    hi = 1.5;
    solver.setBounds(nbd,lo,hi);
#else
    CGfast solver;
    int restartStride = model.dof()/2;
    int printStride = 100;//restartStride;
    double tolLS = 1.0+0*pgtol/sqrt(model.dof());
    int maxIterLS = 20;
    double sigma = 5.0e-6;
    solver.setParameters(CGfast::PR,
			 maxIter,restartStride,printStride,
			 pgtol,pgtol,tolLS,maxIterLS,sigma);

#endif
    // add some viscosity 
    double viscosity = 0.0e-6*scalingFactor;
    ViscousRegularizer vr(bd.nodes(), viscosity);
    bd.pushBack( &vr ); 
  
    double vrTol = 1.0e-10;


//     // deactivate elements 
//   ////////////////////////////////////////////////////////////////////
//   // Skip elements on some of the dodecahedral faces.
//   // 
//   // 1. Check to make sure that the mesh appears to be of a
//   // dodecahedron.  Number of elements should be an integer multiple
//   // of 12, i.e., nElements=n*12.
//   //
//   float triPerFace=0.0;
//   float remainder=0.0;
//   remainder = modf( ntri/12.0, &triPerFace );
//   if( remainder !=0.0 ) {
//     std::cout << "ntri/12 = " << ntri/12.0 << " is not an integer."
// 	      << std::endl;
//     return 0;
//   }
//   // 
//   // 2. Deactivate the first n=nElements/12 elements.
//   //
//   vector<int> removed;
//   // to remove face N add below
//   //      removed.push_back(N);
//   removed.push_back(0);

//   // loop through pentamers to be removed
//   for(int r=0; r<removed.size(); r++) {
//     int face = removed[r];
//     // loop through triangles on each pentamer
//     for( int e=0; e<triPerFace; e++ ) {
//       int tri = triPerFace*face + e;
//       // deactivate triangle in each body
//       for(int b=0; b<bdc.size(); b++) 
// 	bdc[b]->deactivate( tri );
//       //
//       // if springs have been added, deactivate those too, 3 per triangle
//       //
//       if(springsON) {
// 	  bd.deactivate(ntri + 3*tri + 0);
// 	  bd.deactivate(ntri + 3*tri + 1);
// 	  bd.deactivate(ntri + 3*tri + 2);
//       }
//     }      
//   }
  

  for(int n=0; n<nodes.size(); n++) {
    for(int i=0; i<nodes[n]->dof(); i++) nodes[n]->setForce(i,0.0);
  }
    
  for(int b=0; b<bdc.size(); b++) {
    bdc[b]->compute(true,true,false);    
  }
  model.print("deactivated");

  //  
  ////////////////////////////////////////////////////////////////////

  // print starting energy, minimize, and output relaxed shape
  std::cout << "Relaxing shape for gamma = " << gamma << std::endl
	    << "Energy = " << solver.function()/scalingFactor << std::endl;    

    int viter=0;
    double vrEnergy = vr.energy();
    double totalEnergy = solver.function();
    do {
      if(verbose) std::cout << std::endl 
			    << "VISCOUS ITERATION: " << viter 
			    << "\t viscosity = " << vr.viscosity()
			    << std::endl
			    << std::endl;
      viter++;
      vr.step();
      solver.solve( &model );
      double bendingEnergy = bd.energy();
      double stretchEnergy = 0.0;
      if(bdc.size()==2) stretchEnergy = bdc[1]->energy();
      vrEnergy = vr.energy();
      totalEnergy = solver.function();

      if(verbose) {
	std::cout << "viscous energy = " << vrEnergy << std::endl
		  << "bending energy = " << bendingEnergy << std::endl
		  << "stretch energy = " << stretchEnergy << std::endl
		  << "  total energy = " << solver.function() << std::endl
		  << std::endl;
      }
    } while(vrEnergy > vrTol * totalEnergy);
  
  
//     for(int b=0; b<bdc.size(); b++) {
//       Model tempModel(nodes);
//       tempModel.pushBackBody( bdc[b] );
//       std::cout << std::endl;
//       std::cout << "Consistency for body " <<  b << std::endl;
//       std::cout << std::endl;
//       tempModel.checkConsistency(true,false);
//     }
  
    std::cout << "Done relaxing partial capsid." 
	      << std::endl << std::endl
	      << "\t Energy = " << solver.function()/scalingFactor 
	      << std::endl
	      << std::endl;
    char fname[100];
    sprintf(fname,"%s-%g",modelName.c_str(), gamma);
    model.print(fname);
    

  std::cout << "All done." << std::endl;
  return 0;
}
