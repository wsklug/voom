#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include <tvmet/Vector.h>
#include "Node.h"
#include "GLElastic.h"
#include "C0MembraneGL.h"
#include "MembraneGLImplicitMass.h"
#include "GenericBody.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "TriangleQuadrature.h"
#include "C0MembraneGLoutput.h"
#include "S2.h"

#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkWarpVector.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkGeometryFilter.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

struct InputFile {

  // Constructor implemented below...
  InputFile(const string & fname);
  
  double R;
  double r;
  double Y;
  double nu;
  double g0;
  double dg;
  double Gamma;
  double muA;
  double Mx0;
  double Mx1;
  double Mrho0;
  double dt;
  double bfgs_factr;
  double bfgs_pgtol;
  
  int steps;
  int bfgs_iprint;
  int bfgs_m;
  int bfgs_maxIter;

  string vtk_file;
  string topology;

  bool verbose;
};

int main(int argc, char* argv[])
{
  if( argc != 2 ) {
      cout << "Usage: grow modelName." << endl
	   << endl << endl;
	   
      return(0);
  }

  // Alert user to OpenMP status
#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  ////////////////////////////////////////////////////////////////////
  //
  // Input section - load in a VTK file, create nodes, map positions
  // to sphere, and copy element connectivities
  //
  ////////////////////////////////////////////////////////////////////

  string modelName = argv[1];

  // load input file
  const InputFile parameters(modelName);

  ////////////////////////////////////////////////////////////////////
  // read in vtk file
  ////////////////////////////////////////////////////////////////////

//   vtkDataSetReader * reader = vtkDataSetReader::New();
//   reader->SetFileName( parameters.vtk_file.c_str() );
//   // send through normals filter to ensure that triangle orientations
//   // are consistent
//   vtkGeometryFilter * poly = vtkGeometryFilter::New();
//   poly->SetInput( reader->GetOutput() );
//   vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
//   normals->SetInput( poly->GetOutput() );
//   normals->ConsistencyOn();
//   normals->SplittingOff();
//   normals->AutoOrientNormalsOn();
//   vtkPolyData * mesh = normals->GetOutput();
  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( parameters.vtk_file.c_str() );
  vtkDataSet * mesh = reader->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;
  std::cout << "mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  ////////////////////////////////////////////////////////////////////
  // create vector of nodes, get positions from vtk input
  ////////////////////////////////////////////////////////////////////
  int dof=0;
  Model::NodeContainer nodes;
  C0MembraneGL::DefNodeContainer defNodes;
  C0MembraneGL::GLNodeContainer glNodes;

  double Ravg = 0;

  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    Vector3D x;
    mesh->GetPoint(a, &(x[0]));
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(a,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  if(parameters.verbose) cout << "Number of defNodes: " <<defNodes.size() << endl
		   << "Ravg = " << Ravg << endl;

  ////////////////////////////////////////////////////////////////////
  // Map nodes on to sphere.  Two options here: 
  //
  // 1. Input mesh is already a hemisphere.  In this case, treat the
  //    input mesh as the current/deformed configuration, and use
  //    equal-area spherical projection to map the nodes onto a flat
  //    disc representing the stress-free reference configuration.
  //
  // 2. Input mesh has the topology of a complete sphere (maybe an
  //    icosahedron). In this case, treat the input mesh as the
  //    reference configuration and project the nodes radially onto a
  //    sphere to generate the deformed configuration.
  //
  // Note: points are scaled as necessary to ensure that the sphere
  //       has a radius, R.
  ////////////////////////////////////////////////////////////////////

  if( parameters.topology == "hemisphere" ) {
    // Input mesh is of a hemisphere of radius 1, with apex at z=1
    //
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->position();

      x*= parameters.R/Ravg;

      defNodes[i]->setPoint(x);

      // Compute reference configuration by equal-area spherical projection
      
      // apex
      Vector3D x0;
      x0 = 0.0, 0.0, parameters.R;
      
      double r = norm2(x-x0);
      x(2) = 0.0;
      x = (r/norm2(x))*x;
      x(2) = parameters.R;
      defNodes[i]->setPosition(x);

    }

  } else if (parameters.topology == "sphere") {

    // Assume input mesh is of spherical topology (maybe an icosahedron)
    
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->position();
      x *= parameters.R/Ravg;
      defNodes[i]->setPosition(x);
      
      // initially spherical
      //x = defNodes[i]->point();
      x *= parameters.R/tvmet::norm2(x);
      defNodes[i]->setPoint(x);
    }
  } else {
    cout << "\"" << parameters.topology 
	 << "\" is an unknown topology type.  Exiting." 
	 << endl;
    return 0;
  }

  ////////////////////////////////////////////////////////////////////
  // Create and initialize concentration nodes.  If within a radius of
  // r, set concentration to 1.  
  ////////////////////////////////////////////////////////////////////
  for(int a=0; a<defNodes.size(); a++) {
    NodeBase::DofIndexMap idx(1);
    idx[0]=dof++;
    ScalarFieldNode<3>* s = 
      new ScalarFieldNode<3>(a, idx, defNodes[a]->position(), 0.0);
    nodes.push_back( s );
    glNodes.push_back( s );
  }

  for(int a=0; a<defNodes.size(); a++) {
    const Vector3D & x = defNodes[a]->point();
    double r = sqrt(x(0)*x(0)+x(1)*x(1));
    if(x(2) > 0.0 && r < parameters.r ) {
      glNodes[a]->setPoint(0,1.0);
    } else {
      glNodes[a]->setPoint(0,0.0);
    }

  }

  if(parameters.verbose) cout << "Number of glNodes: " <<glNodes.size() << endl
		   << "Number of nodes: " <<nodes.size() << endl;


  ////////////////////////////////////////////////////////////////////
  // Solution section
  ////////////////////////////////////////////////////////////////////

  // initialize material response

  cout << "Creating GLElastic material." << endl;

  double KC = 0.0, KG = 0.0, C0 = 0.0;
    
  int formulation=0;

  GLElastic mat( KC, KG, C0, 
		 parameters.Y, parameters.nu, 
		 parameters.Gamma, parameters.g0, 
		 parameters.dg, parameters.muA, 
		 formulation );

  // create Body 
  GenericBody bd;

  cout << "Creating ShapeTri3." << endl;


  int ntri=mesh->GetNumberOfCells();
  if(parameters.verbose) cout << "Number of triangles: " <<ntri << endl;
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

  TriangleQuadrature quad(quadOrder);
  TriangleQuadrature quad2(quadOrder+1);


  cout << "Creating elements." << endl;

  std::vector< C0MembraneGL* > membranes;
  std::vector< MembraneGLImplicitMass* > masses;

  ////////////////////////////////////////////////////////////////////
  // copy connectivities
  ////////////////////////////////////////////////////////////////////
  for (int e = 0; e<ntri; e++){

    C0MembraneGL::DefNodeContainer dnodes(nodes_per_element);
    C0MembraneGL::GLNodeContainer gnodes(nodes_per_element);
    
    for(int a=0; a<nodes_per_element; a++) {
      int A = mesh->GetCell(e)->GetPointId(a);
      dnodes[a] = defNodes[A];
      gnodes[a] = glNodes[A];
    }
    
    C0MembraneGL* elem = new C0MembraneGL(dnodes, gnodes, &mat, &quad, shape);
    membranes.push_back( elem );
    bd.addElement( elem );
    
    if( nodes_per_element == 3 ) {
      
      MembraneGLImplicitMass * mass 
	= new MembraneGLImplicitMass(dnodes, gnodes, &quad2, shape, 
				     parameters.Mx0, parameters.Mx1, 
				     parameters.Mrho0);
      masses.push_back( mass );
      bd.addElement( mass );

    } else if( nodes_per_element == 6 ) {
      // break element up into 4 3-noded elements to compute masses

      int idx[4][3] = { 0,3,5, 
			1,4,3, 
			2,5,4, 
			3,4,5 };
      for(int t=0; t<4; t++) {
	dnodes.resize(3);
	gnodes.resize(3);
	
	for(int a=0; a<3; a++) {
	  int A = mesh->GetCell(e)->GetPointId(idx[t][a]);
	  dnodes[a] = defNodes[A];
	  gnodes[a] = glNodes[A];
	}
	
	MembraneGLImplicitMass * mass 
	  = new MembraneGLImplicitMass(dnodes, gnodes, &quad, shape3, 
				       parameters.Mx0, parameters.Mx1, 
				       parameters.Mrho0);
	masses.push_back( mass );
	bd.addElement( mass );
	
      }
    }
  }

  for( int a=0; a<nodes.size(); a++ ) bd.addNode( nodes[a] );

  for(int a=0; a<nodes.size(); a++) 
    for(int i=0; i<nodes[i]->dof(); i++) nodes[i]->setForce(i,0.0);

  bd.compute(true,true,false);

  cout << "Creating Model." << endl;

  Model model(nodes);

  model.pushBackBody( &bd );

  // initialize Quasi-Newton BFGS solver
  Lbfgsb solver(model.dof(), parameters.bfgs_m, 
		parameters.bfgs_factr, parameters.bfgs_pgtol, 
		parameters.bfgs_iprint, parameters.bfgs_maxIter );

  // set up bounds for solver
  bool boundaryConditionFlag = true;
  blitz::Array<int,1> nbd(model.dof());
  blitz::Array<double,1> lo(model.dof());
  blitz::Array<double,1> hi(model.dof());
  nbd = 0;
  lo = 0;
  hi = 0;
  
  if(parameters.Y <= 0.0 ) {
    // fix all deformation nodes
    cout << endl << "Fixing all deformation nodes." << endl << endl;
    for(int a=0, ia=0; a<defNodes.size(); a++) {
      for(int i=0; i<3; i++, ia++) {
	nbd(ia) = 2;
	hi(ia) = lo(ia) = defNodes[a]->getPoint(i);
      }
    }
  }
  
  // bound all gl nodes
  for(int ia=3*defNodes.size(); ia<3*defNodes.size()+glNodes.size(); ia++) {
    nbd(ia) = 2;
    lo(ia) = 0.0;
    hi(ia) = 1.0;
  }
  
  solver.setBounds(nbd,lo,hi);
  
  // constrain all nodes to lie on the sphere
  for(int a=0; a<defNodes.size(); a++) {
    S2 * s = new S2( defNodes[a], parameters.R );
    bd.addConstraint( s );
  }
  
  for(int m=0; m<masses.size(); m++) {
    masses[m]->step(parameters.dt);
  }

  for(int n=0; n<nodes.size(); n++) {
    for(int i=0; i<nodes[n]->dof(); i++) nodes[n]->setForce(i,0.0);
  }
    
//   bd.compute(true,true,false);    
  model.computeAndAssemble(solver, true, true, false);

  //  
  ////////////////////////////////////////////////////////////////////

  C0MembraneGLoutput output;
  char fname[100];
  sprintf(fname,"%s-0",modelName.c_str());
  output(membranes, defNodes, glNodes, fname);
  
  int step=1;
  do {
    if(parameters.verbose) std::cout << std::endl 
			  << "time step: " << step
			  << std::endl
			  << std::endl;
    for(int m=0; m<masses.size(); m++) {
      masses[m]->step(parameters.dt);
    }
    solver.solve( &model );
    
    if(parameters.verbose) {
      std::cout << "  total energy = " << solver.function() << std::endl
		<< std::endl;
    }
    
    
    sprintf(fname,"%s-%d",modelName.c_str(), step);
    output(membranes, defNodes, glNodes, fname);
    
//     bd.checkConsistency();

    step++;
    
  } while(step < parameters.steps);//vrEnergy > vrTol * totalEnergy);
  
  std::cout << "All done." << std::endl;
  return 0;
}

InputFile::InputFile(const string & fname) {

  ifstream ifs( fname.c_str() );
  if (!ifs) {
    cout << endl
	 << "Cannot open file " << fname
	 << endl
	 << endl;
    exit (0);
  }

  R = 0.0;
  r = 0.0;
  Y = 0.0;;
  nu = 0.0;
  g0 = 0.0;
  dg = 0.0;
  Gamma = 0.0;
  muA = 0.0;
  Mx0 = 0.0;
  Mx1 = 0.0;
  Mrho0 = 0.0;
  dt = 0.0;
  bfgs_factr = 0.0;
  bfgs_pgtol = 0.0;
  
  bfgs_iprint = 0;
  steps = 0;
  bfgs_m = 0;
  bfgs_maxIter = 0;

  vtk_file = "";
  topology = "";

  // read through file line by line and assign parameter values
  string line;
  int number;
  string name;
  cout << "Parsing input file.  List of (key, value) pairs:" << endl;
  while ( getline( ifs, line ) ) {
    string::size_type pos=line.find(":");
    if( pos == string::npos ) continue;
    
    string key = line.substr(0,pos);
    string value = line.substr(pos+1);

    // strip white space from value
    for(pos=value.find(" "); pos != string::npos; pos=value.find(" ") ) {
      value.erase(pos,1);
    }
    for(pos=value.find("\t"); pos != string::npos; pos=value.find("\t") ) {
      value.erase(pos,1);
    }

    cout << "\"" << key << "\", \""<< value << "\""<< endl;

    // assign values 
    if( key == "vtk file" ) {	
      vtk_file = value;
    } else if( key == "sphere radius" ) {
      R = atof(value.c_str());
    } else if( key == "nucleation radius" ) {
      r = atof(value.c_str());
    } else if( key == "topology" ) {
      topology = value;
    } else if( key == "Young\'s modulus" ) {
      Y = atof(value.c_str());
    } else     if( key == "Poisson's ratio" ) {
      nu = atof(value.c_str());
    } else if( key == "GL curvature" ) {
      g0 = atof(value.c_str());
    } else if( key == "GL deltaG" ) {	
      dg = atof(value.c_str());
    } else if( key == "Gamma" ) {
      Gamma = atof(value.c_str());
    } else if( key == "muA" ) {
      muA = atof(value.c_str());
    } else if( key == "Mx0" ) {
      Mx0 = atof(value.c_str());
    } else if( key == "Mx1" ) {
      Mx1 = atof(value.c_str());
    } else if( key == "Mrho0" ) {		
      Mrho0 = atof(value.c_str());
    } else if( key == "time step" ) {	
      dt = atof(value.c_str());
    } else if( key == "number steps" ) {	
      steps = atoi(value.c_str());
    } else if( key == "BFGS iprint" ) {	
      bfgs_iprint = atoi(value.c_str());
    } else if( key == "BFGS factr" ) {	
      bfgs_factr = atof(value.c_str());
    } else if( key == "BFGS pgtol" ) {	
      bfgs_pgtol = atof(value.c_str());
    } else if( key == "BFGS m" ) {	
      bfgs_m = atoi(value.c_str());
    } else if( key == "BFGS maxIter" ) {
      bfgs_maxIter = atoi(value.c_str());
    } else if( key == "verbose" ) {
      verbose = atoi(value.c_str());
    } else { 
      cout << "\"" << key << "\"" 
	   << " is an undefined keyword.  Did you spell it incorrectly?"
	   << endl;
      exit(0);
    }
  }
  return;
}
