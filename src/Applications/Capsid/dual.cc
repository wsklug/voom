#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"

#include <vtkCell.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void computeAspherity(char * modelName, double & aspherity, double & radius);

int main(int argc, char* argv[])
{
  if( argc != 2 ) {    
    cout << "argc = " << argc << endl
	 << "Usage: dual modelName" << endl;
      return(0);
  }

  ////////////////////////////////////////////////////////////////////
  // Input section
  ////////////////////////////////////////////////////////////////////

  bool debug=true;

  // read in vtk file
  string inputFileName = argv[1];  
  int found = inputFileName.find(".vtk");

  string modelName;
  if(found!=std::string::npos){ 
    modelName = inputFileName.substr(0,found);
  }

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  //We will use this object, shortly, to ensure consistent triangle orientations
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();

  //We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput()
  //If our input vtk file has vtkUnstructuredGridData instead of vtkPolyData
  //then we need to convert it using vtkGeometryFilter
  vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();
  ds->Update();
  if(ds->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = 
      reader->GetUnstructuredGridOutput();    
    vtkSmartPointer<vtkGeometryFilter> geometryFilter = 
      vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInput(unstructuredGrid);
    geometryFilter->Update(); 
    vtkSmartPointer<vtkPolyData> polydata = geometryFilter->GetOutput();
    normals->SetInput( polydata);
  }
  else{
    normals->SetInput(reader->GetOutput());
  }

  // send through normals filter to ensure that triangle orientations
  // are consistent
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
//   normals->Update();
  vtkPolyData * vtkMeshOld = normals->GetOutput();
  vtkMeshOld->Update();
  int npts = vtkMeshOld->GetNumberOfPoints(); 
  std::cout << "npts = " 
	    << npts
	    << std::endl;

  // get vertex positions
  std::vector< tvmet::Vector<double,3> > points( npts );
  for(int a=0; a<npts; a++) {
    vtkMeshOld->GetPoint(a, &(points[a](0)));
  }

  //Can we add displacements to vertex positions before generatinng
  //the dual mesh - Amit
  
  //get displacements
  //std::vector< tvmet::Vector<double,3> > displacements( npts );
  string vectorName="displacements";
  vtkSmartPointer<vtkDataArray> displacements = vtkMeshOld->GetPointData()->
    GetVectors(vectorName.c_str());
  double currDisp[3];
  for(int a=0; a<npts; a++){
    displacements->GetTuple(a,currDisp);
    points[a][0] += currDisp[0];
    points[a][1] += currDisp[1];
    points[a][2] += currDisp[2];
  }


  // get connectivities
  int ntriOld=vtkMeshOld->GetNumberOfCells();  
  HalfEdgeMesh::ConnectivityContainer connectivitiesOld( ntriOld );
  cout << "Number of triangles: " <<ntriOld << endl;
  for (int i = 0; i<ntriOld; i++){
    assert(vtkMeshOld->GetCell(i)->GetNumberOfPoints() == 3);
    for(int a=0; a<3; a++) {
      connectivitiesOld[i](a) = vtkMeshOld->GetCell(i)->GetPointId(a);
    }
  }

  // Create new points for refined mesh
  //
  // new point at barycenter of each face
  std::vector< tvmet::Vector<double,3> > dualPoints( connectivitiesOld.size() );
  for(int f=0; f<connectivitiesOld.size(); f++) {
    tvmet::Vector<double,3> x(0.0);
    for(int v=0; v<3; v++) {
      x += points[connectivitiesOld[f](v)];
    }
    x/=3;
    dualPoints[f] = x;
  }

  
  HalfEdgeMesh primalMesh(connectivitiesOld,npts);



  // build new connectivities, one face for each primal vertex
  int nDualFaces = points.size();
  std::vector< std::vector< int > > dualConnectivities(nDualFaces);


  if(debug) {
    cout << endl << "Building connectivities for " << nDualFaces << " new faces."
	 << endl;
  }

  int nFaceVertices = 0;
  for(int f=0; f<nDualFaces; f++) {
    
    if(debug) cout << endl << "Face " << f << endl;

    dualConnectivities[f].resize( primalMesh.vertices[f]->halfEdges.size() );
    // walk around this primal vertex's halfedges in CCW order, and
    // get the id of each one's face
    HalfEdge* H=primalMesh.vertices[f]->halfEdges[0];
    for(int v=0; v<dualConnectivities[f].size(); v++) {
      
      dualConnectivities[f][v] = H->face->id;
      H = H->opposite->prev;

      nFaceVertices++;

      if(debug) {
	cout << dualConnectivities[f][v] << " ";
      }
    }
    if(debug) cout << endl;

  }   

  // new triangulate the dual mesh
  std::vector< tvmet::Vector<double,3> > dualTriPoints( dualPoints );
  
  // to copy old vertices to centers of dual faces
//   for(int v=0; v<points.size(); v++) {
//     dualTriPoints.push_back( points[v] );
//   }

  // to compute barycenters of dual faces
  for(int f=0; f<dualConnectivities.size(); f++) {
    tvmet::Vector<double,3> x(0.0);
    for(int v=0; v<dualConnectivities[f].size(); v++) {
      x += dualPoints[dualConnectivities[f][v]];
    }
    x/=dualConnectivities[f].size();
    dualTriPoints.push_back( x );
  }

  std::vector<int> valences;
  std::vector< std::vector< int > > dualTriangles;
  for(int f=0; f<dualConnectivities.size(); f++) {
    int nvf=dualConnectivities[f].size();
    for(int v=0; v<nvf; v++) {
      std::vector<int> vf(3);
      vf[0]=dualConnectivities[f][v];
      vf[1]=dualConnectivities[f][(v+1)%nvf];
      vf[2]=dualPoints.size()+f ; // center of the face
      dualTriangles.push_back(vf);
      valences.push_back(dualConnectivities[f].size());
    }
  }
  

  // print out new mesh to a vtk file
  char  fileName[50];
  sprintf(fileName, "%s-dual.vtk",modelName.c_str());
  cout << "Saving dual mesh as " << fileName << endl;

  std::ofstream ofs(fileName);
  if (!ofs) {
    std::cout << "can not open output ("
	      << fileName
	      << ") file." << std::endl;
    exit(0);
  }
  
  // print vertices
    ofs << "# vtk DataFile Version 2.0\n"
	<< modelName << " dual "
	<< std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << dualPoints.size() << "  double" << std::endl;
    ofs.setf(std::ios_base::scientific,std::ios_base::floatfield);
    for (int i=0; i<dualPoints.size(); i++) {
      ofs << std::setprecision(16) 
	  << std::setw(24) << dualPoints[i](0) 
	  << std::setw(24) << dualPoints[i](1) 
	  << std::setw(24) << dualPoints[i](2) 
	  << std::endl;
    }

    // print faces
    ofs << "POLYGONS  " << dualConnectivities.size() << "  "
	<< dualConnectivities.size() + nFaceVertices << std::endl;
    for (int f=0; f<dualConnectivities.size(); f++) {
      ofs << dualConnectivities[f].size() << "  ";
      for(int v=0; v<dualConnectivities[f].size(); v++) {
	ofs << std::setw(10) << dualConnectivities[f][v];
      }
      ofs << std::endl;
    }


    ofs << "CELL_DATA    " << dualConnectivities.size() << std::endl;
    //
    // output for strain energy
    ofs << "SCALARS    valence    int    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for (int f=0; f<dualConnectivities.size(); f++) {
      ofs << dualConnectivities[f].size() << std::endl;
    }
    ofs << std::endl;



  sprintf(fileName, "%s-dual-tri.vtk",modelName.c_str());
  cout << "Saving triangulated dual mesh as " << fileName << endl;

  std::ofstream ofs_tri(fileName);
  if (!ofs_tri) {
    std::cout << "can not open output ("
	      << fileName
	      << ") file." << std::endl;
    exit(0);
  }
  
  // print vertices
    ofs_tri << "# vtk DataFile Version 2.0\n"
	<< modelName << " dual "
	<< std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << dualTriPoints.size() << "  double" << std::endl;
    ofs_tri.setf(std::ios_base::scientific,std::ios_base::floatfield);
    for (int i=0; i<dualTriPoints.size(); i++) {
      ofs_tri << std::setprecision(16) 
	  << std::setw(24) << dualTriPoints[i](0) 
	  << std::setw(24) << dualTriPoints[i](1) 
	  << std::setw(24) << dualTriPoints[i](2) 
	  << std::endl;
    }

    // print faces
    ofs_tri << "POLYGONS  " << dualTriangles.size() << "  "
	    << 4*dualTriangles.size() << std::endl;
    for (int f=0; f<dualTriangles.size(); f++) {
      ofs_tri << dualTriangles[f].size() << "  ";
      for(int v=0; v<dualTriangles[f].size(); v++) {
	ofs_tri << std::setw(10) << dualTriangles[f][v];
      }
      ofs_tri << std::endl;
    }

    ofs_tri << "CELL_DATA    " << dualTriangles.size() << std::endl;
    //
    // output for strain energy
    ofs_tri << "SCALARS    valence    int    1" << std::endl;
    ofs_tri << "LOOKUP_TABLE default" << std::endl;
    for (int f=0; f<dualTriangles.size(); f++) {
      ofs_tri << valences[f] << std::endl;
    }
    ofs_tri << std::endl;

    return 0;
}

