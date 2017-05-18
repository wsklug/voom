#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>
#include <vtkWindowedSincPolyDataFilter.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;

int main(int argc, char* argv[])
{
  if( argc != 3 ) {
      cout << "Usage: indent modelName decimationTarget" << endl;
      return(0);
  }

  string modelName = argv[1];

   string inputFileName = modelName + ".vtk";

   // read in file
   vtkDataSetReader * reader = vtkDataSetReader::New();
   reader->SetFileName( inputFileName.c_str() );

//  string inputFileName = modelName + ".vtu";

//  // read in file
//  vtkUnstructuredGridReader * reader = vtkUnstructuredGridReader::New();
//  reader->SetFileName( inputFileName.c_str() );

  // extract boundary surface
  vtkGeometryFilter * geometry = vtkGeometryFilter::New();
  geometry->SetInput( reader->GetOutput() );

  // send through normals filter to ensure that triangle orientations
  // are consistent
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  normals->SetInput( geometry->GetOutput() );
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();

  // get the mesh
  vtkPolyData * mesh = normals->GetOutput();
  mesh->Update();
  
  std::cout << "Surface Mesh: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // convert to triangle mesh if necessary
  vtkTriangleFilter * triangles = vtkTriangleFilter::New();
  triangles->SetInput( normals->GetOutput() );
  triangles->Update();

  mesh = triangles->GetOutput();
  mesh->Update();
  std::cout << "Surface Triangulation: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // write basic triangle mesh to a file
  vtkPolyDataWriter * writer = vtkPolyDataWriter::New();
  writer->SetInput( triangles->GetOutput() );
  
  string outputFileName = modelName + "-surface.vtk";
  writer->SetFileName( outputFileName.c_str() );
  writer->SetFileTypeToBinary();
//   writer->SetFileTypeToASCII();
  writer->Write();

  // smooth using vtkWindowedSincPolyDataFilter
  vtkWindowedSincPolyDataFilter * smoother=vtkWindowedSincPolyDataFilter::New();
  smoother->SetInput( triangles->GetOutput() );
  smoother->SetNumberOfIterations( 20 );
  smoother->BoundarySmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->SetFeatureAngle( 150.0 );
  smoother->SetPassBand( 0.001 );

  outputFileName = modelName + "-smoothed.vtk";
  writer->SetInput( smoother->GetOutput() );
  writer->SetFileName( outputFileName.c_str() );
  writer->Write();

  // decimate using vtkDecimatePro
  //
  //  http://www.vtk.org/doc/release/5.4/html/a00397.html
  //
  // Note: two other similar decimation filters are 
  //  vtkQuadricClustering: http://www.vtk.org/doc/release/5.4/html/a01444.html
  //  vtkQuadricDecimation: http://www.vtk.org/doc/release/5.4/html/a01446.html
  //
  vtkDecimatePro * decimate = vtkDecimatePro::New();
//   decimate->SetInput( triangles->GetOutput() );
  decimate->SetInput( smoother->GetOutput() );
  // How much to decimate, i.e., 0.9 tries to remove 90% of the triangles
  double target = atof(argv[2]);
  decimate->SetTargetReduction( target );
  decimate->PreserveTopologyOn();
  decimate->SetFeatureAngle(179.0);
  decimate->SplittingOff();
  // check vtk online documentation for many many other options for
  // this filter
  decimate->Update();
  
  mesh = decimate->GetOutput();
  mesh->Update();
  std::cout << target*100.0 << "\% Decimated Triangulation: " << std::endl
	    << " mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << " mesh->GetNumberOfCells() = " << mesh->GetNumberOfCells()
	    << std::endl;

  // write decimated mesh to file
  writer->SetInput( decimate->GetOutput() );
  char tmp[20];
  sprintf(tmp,"%.2g",target);
  outputFileName = modelName + "-decimated-" + tmp + ".vtk";
  writer->SetFileName( outputFileName.c_str() );
  writer->Write();

  // smooth decimated mesh
  vtkWindowedSincPolyDataFilter * second_smoother
    =vtkWindowedSincPolyDataFilter::New();
  second_smoother->SetInput( decimate->GetOutput() );
  second_smoother->SetNumberOfIterations( 20 );
  second_smoother->BoundarySmoothingOn();
  second_smoother->NormalizeCoordinatesOn();
  second_smoother->SetFeatureAngle( 150.0 );
  second_smoother->SetPassBand( 0.001 );
  second_smoother->SetInput( decimate->GetOutput() );

  // write smoothed decimated mesh to file
  writer->SetInput( second_smoother->GetOutput() );
  outputFileName = modelName + "-decimated-" + tmp + "-smoothed.vtk";
  writer->SetFileName( outputFileName.c_str() );
  writer->Write();

  mesh = second_smoother->GetOutput();
  mesh->Update();

  // write nodes and connectivities of decimated mesh to files
  outputFileName = modelName + ".nd";
  ofstream node_file(outputFileName.c_str() );
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    for(int i=0; i<3; i++) {
      node_file << mesh->GetPoint(a)[i] << " ";
    }
    node_file << std::endl;
  }

  outputFileName = modelName + ".el";
  ofstream element_file(outputFileName.c_str() );
  for(int e=0; e<mesh->GetNumberOfCells(); e++) {
    for(int a=0; a<3; a++) {
      element_file << mesh->GetCell(e)->GetPointId(a) << " ";
    }
    element_file << std::endl;
  }

  return 0;
}

