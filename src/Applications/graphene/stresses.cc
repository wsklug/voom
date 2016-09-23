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
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

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

  int nSteps = std::atoi(argv[2]);

  string inputFileName = modelName + "-step0000.vtk";
  vtkSmartPointer<vtkDataSetReader> reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  // // send through normals filter to ensure that triangle orientations
  // // are consistent
  // vtkSmartPointer<vtkPolyDataNormals> normals = vtkPolyDataNormals::New();
  // normals->SetInput( reader->GetOutput() );
  // normals->ConsistencyOn();
  // normals->SplittingOff();
  // normals->AutoOrientNormalsOn();
  vtkSmartPointer<vtkDataSet> mesh = reader->GetOutput();
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

  // // read in triangle connectivities
  // vector< tvmet::Vector<int,3> > connectivities;
  // tvmet::Vector<int, 3> c;
  // int ntri=mesh->GetNumberOfCells();
  // connectivities.reserve(ntri);
  // if(verbose) cout << "Number of triangles: " <<ntri << endl;
  // for (int i = 0; i<ntri; i++){
  //   assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
  //   for(int a=0; a<3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
  //   connectivities.push_back(c);
  // }

  	
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
  Graphene stretching;

  if( nodes_per_element == 3 ) {
    typedef C0MembraneBody<TriangleQuadrature,Graphene,ShapeTri3> MB;
    MB * bdm = new MB(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm->setOutput(paraview);
  
    // Vector3D N;
    // N = 1.0, 0.0, 0.0;
    // bdm->addStressDirection(N);
    // N = 0.0, 1.0, 0.0;
    // bdm->addStressDirection(N);
    
    bdc.push_back(bdm);

  } else {
    typedef C0MembraneBody<TriangleQuadrature,Graphene,ShapeTri6> MB;
    MB * bdm = new MB(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm->setOutput(paraview);
  
    // Vector3D N;
    // N = 1.0, 0.0, 0.0;
    // bdm->addStressDirection(N);
    // N = 0.0, 1.0, 0.0;
    // bdm->addStressDirection(N);
    
    bdc.push_back(bdm);
  }

  
  string zfsName = modelName + ".zfs";
  ofstream zfs(zfsName.c_str());

  std::cout << setw(6) << "Step" 
	    << setw(12) << "zmax"
	    << setw(12) << "load"
	    << setw(12) << "restraint"
	    << setw(12) << "s0min"
	    << setw(12) << "s0max"
	    << setw(12) << "s1min"
	    << setw(12) << "s1max" 
	    << std::endl;

  for(int step=0; step<=nSteps; step++) {

    // std::cout << "Step " << step << std::endl;

    char name[100]; 
    sprintf(name,"%s-step%04d.vtk",modelName.c_str(),step);

    vtkSmartPointer<vtkDataSetReader> step_reader = vtkDataSetReader::New();
    step_reader->SetFileName( name );
    vtkSmartPointer<vtkDataSet> step_mesh = step_reader->GetOutput();
    
    // std::cout << "step_mesh created." << std::endl;

    step_mesh->Update();

    // std::cout << "Mesh from " << name << " updated." << std::endl;

    // std::cout << "Point data contains " 
    // 	      << step_mesh->GetPointData()->GetNumberOfArrays()
    // 	      << " arrays:" << std::endl;

    // for(int d=0; d<step_mesh->GetPointData()->GetNumberOfArrays(); d++) 
    //   std::cout << step_mesh->GetPointData()->GetArrayName(d) << " ";
    // std::cout << std::endl;

    vtkDataArray * displacements 
      = step_mesh->GetPointData()->GetVectors("displacements");

    double zmax = std::abs( displacements->GetComponent(0,2) );

    for(int a=0; a<nodes.size(); a++) {
      // initialize nodal forces
      for(int i=0; i<nodes[a]->dof(); i++) nodes[a]->setForce(i,0.0);
      
      // apply displacements from vtk file
      Vector3D x(0.0);
      displacements->GetTuple( a, &(x(0)) );
      
      x += defNodes[a]->position();
      defNodes[a]->setPoint(x);
      
      // look for max z
      zmax = std::max( zmax, std::abs(  displacements->GetComponent(0,2) ) );
    }
    
    // std::cout << "zmax = " << zmax << std::endl;
    
    // step_mesh->Delete();
    // step_normals->Delete();
    // step_reader->Delete();
    
    // Compute stresses etc.
    bdc[0]->compute(true,true,false);    
    
    // compute out of balance force (indentation force) in z direction
    double load=0.0;
    double restraint=0.0;
    for(int a=0; a<nodes.size(); a++) {
      if( tvmet::norm2( defNodes[a]->position() ) < 0.9*Rmax ) 
	load += nodes[a]->getForce(2); 
      else
	restraint += nodes[a]->getForce(2); 
    }
    
    // std::cout << "residual = " << residual << std::endl;

    // std::cout << "Printing out vtk file with stresses... ";

    // Print stuff out to new vtk file
    sprintf(name,"%s-stress-step%04d",modelName.c_str(),step);
    bdc[0]->printParaview(name);

    // std::cout << "done." << std::endl;


    // Load the new vtk file and find max stresses

    sprintf(name,"%s-stress-step%04d.vtk",modelName.c_str(),step);

    // std::cout << "Reading in vtk file with stresses, "
    // 	      << name << "... ";

    continue;

    vtkSmartPointer<vtkDataSetReader> stress_reader 
      = vtkDataSetReader::New();
    stress_reader->SetFileName( name );
    stress_reader->Update();

    vtkSmartPointer<vtkDataSet> stress_mesh = stress_reader->GetOutput();

    stress_mesh->Update();

    // std::cout << "stress_mesh->GetNumberOfPoints() = " << stress_mesh->GetNumberOfPoints()
    // 	    << std::endl;

    // std::cout << "done." << std::endl;

    // std::cout << "Computing min and max S0... ";
    // std::cout.flush();

    vtkCellData * cd = stress_mesh->GetCellData();

    // std::cout << "Cell data contains " << cd->GetNumberOfArrays()
    // 	      << " arrays:" << std::endl;

    // for(int d=0; d<cd->GetNumberOfArrays(); d++) 
    //   std::cout << cd->GetArrayName(d) << " ";
    // std::cout << std::endl;
    
    stress_mesh->GetCellData()->SetActiveScalars("S0");

    // std::cout << "Set active scalars to S0." << std::endl;

    double s[2];
    stress_mesh->GetCellData()->GetScalars()->GetRange(s);
    double s0min = s[0];
    double s0max = s[1];
    // s0min = s0max = stresses->GetComponent(0,0);
    // for(int e=1; e<ntri; e++) {
    //   s0min = std::min( s0min, stresses->GetComponent(e,0) );
    //   s0max = std::max( s0max, stresses->GetComponent(e,0) );
    // }
    // std::cout << "s0min = " << s0min << std::endl;
    // std::cout << "s0max = " << s0max << std::endl;

    // std::cout << "done." << std::endl;

    // std::cout << "Computing min and max S1... ";

    stress_mesh->GetCellData()->SetActiveScalars("S1");

    // std::cout << "Set active scalars to S1." << std::endl;

    stress_mesh->GetCellData()->GetScalars()->GetRange(s);
    double s1min = s[0];
    double s1max = s[1];

    // std::cout << "s1min = " << s1min << std::endl;
    // std::cout << "s1max = " << s1max << std::endl;

    // std::cout << "done." << std::endl;

    zfs << std::setw( 6 ) << step
	<< std::setw( 24 ) << std::setprecision(16) << zmax
	<< std::setw( 24 ) << std::setprecision(16) << load
	<< std::setw( 24 ) << std::setprecision(16) << restraint
	<< std::setw( 24 ) << std::setprecision(16) << s0min
	<< std::setw( 24 ) << std::setprecision(16) << s0max
	<< std::setw( 24 ) << std::setprecision(16) << s1min
	<< std::setw( 24 ) << std::setprecision(16) << s1max
	<< std::endl;

    std::cout << std::setw( 6 ) << step
	      << std::setw( 12 ) << std::setprecision(4) << zmax
	      << std::setw( 12 ) << std::setprecision(4) << load
	      << std::setw( 12 ) << std::setprecision(4) << restraint
	      << std::setw( 12 ) << std::setprecision(4) << s0min
	      << std::setw( 12 ) << std::setprecision(4) << s0max
	      << std::setw( 12 ) << std::setprecision(4) << s1min
	      << std::setw( 12 ) << std::setprecision(4) << s1max
	      << std::endl;

    // stress_mesh->Delete();
    // stress_normals->Delete();
    // stress_reader->Delete();
    // std::cout << "Deleted vtk objects." << std::endl;
  }
    
  zfs.close();

  std::cout << "Processing complete." << std::endl;

  // mesh->Delete();
  // normals->Delete();
  // reader->Delete();

  return 0;
}

