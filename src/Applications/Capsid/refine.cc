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
#include <vtkDataSet.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void computeAspherity(char * modelName, double & aspherity, double & radius);

int main(int argc, char* argv[])
{
  if( argc != 3 ) {    
    cout << "argc = " << argc << endl
	 << "Usage: refine modelName nSubEdges." << endl;
      return(0);
  }

  ////////////////////////////////////////////////////////////////////
  // Input section
  ////////////////////////////////////////////////////////////////////

  bool debug=false;

  // read in vtk file
  string modelName = argv[1];

  string inputFileName = modelName + ".vtk";

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );

  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();

  //We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput() If
  //our input vtk file has vtkUnstructuredGridData instead of
  //vtkPolyData then we need to convert it using vtkGeometryFilter
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

  HalfEdgeMesh heMeshOld(connectivitiesOld,npts);


  // Create new points for refined mesh
  const int nSub = atoi(argv[2]);
  assert(nSub > 0);
  std::cout << "Refining mesh by subdividing each edge into " 
	    << nSub << " sub-edges." << std::endl;
  
  // 1. Keep points of old mesh.  No action required; just build new
  //    point list from old one.

  // 2. Add new points along edges.  Build a list of the new points on
  //    each half-edge for later when we will build new connectivities.
  int nHE = heMeshOld.halfEdges.size();
  std::vector< std::vector<int> > halfEdgePoints(nHE);
  // visit half-edges by looping through faces
  int nF=heMeshOld.faces.size();
  for(int f=0; f<nF; f++) {   
    const Face * face = heMeshOld.faces[f];
    // for each face loop through half-edges and add vertices if the
    // edge doesn't yet have them
    if(debug) cout << "Splitting halfedges for face " << f << endl;
    for (int h=0; h<3; h++) {
      const HalfEdge * he = face->halfEdges[h];
      int id = he->id;
      if(  halfEdgePoints[id].size() == 0 ) { // edge has no new vertices yet
	// linearly interpolate to get new edge vertex positions
	for(int i=1; i<nSub; i++) {
	  double s=double(i)/nSub;
	  const tvmet::Vector<double, 3> & p1 = points[he->prev->vertex->id];
	  const tvmet::Vector<double, 3> & p2 = points[he->vertex->id];
	  tvmet::Vector<double, 3> p;
	  p = (1.0-s)*p1 + s*p2;
	  halfEdgePoints[id].push_back( npts++ );
	  points.push_back(p);
	}
	// add these new points to the opposite half-edge, in reverse
	// order
	int op = he->opposite->id;
	for(int i=0; i<nSub-1; i++) { 
	  halfEdgePoints[op].push_back( halfEdgePoints[id][nSub-2-i] );
	}
      }
      if (debug) {
	cout << "halfedge " << h << " gets " << halfEdgePoints[h].size() 
	     << " new points" 
	     << endl;
      }
    }
  }
  
  // 3. Add points in the interior of each face, according to the
  // following indexing scheme relying on a multiindex (i0,i1,i2)
  // where i0, etc. are between 0 and nsub.
  //
  // i2
  // :
  // 9   
  // | \   
  // 5 - 8
  // | \ | \       j(i0,i1,i2) = i2+(i1+i2)*(i1+i2+1)/2
  // 2 - 4 - 7   
  // | \ | \ | \      
  // 0 - 1 - 3 - 6... i1    
  //       
  //   i0=0 is halfEdges[2]
  //   i1=0 is halfEdges[0]
  //   i2=0 is halfEdges[1]

  int ptsPerFace=((nSub+1)*(nSub+2))/2;
  std::vector< std::vector<int> > facePoints(nF);
  if(debug) {
    cout << "Each face should get " << ptsPerFace << " total points"
	 << endl;
  }
  for(int f=0; f<nF; f++) {
    const Face * face = heMeshOld.faces[f];
    if(debug) {
      cout << "Adding points to face " << f
	   << endl;
    }
    int i0=nSub, i1=0, i2=0;
    while( i2 <= nSub ) {
      int j = i2+((i1+i2)*(i1+i2+1))/2;
      if(debug) {
	cout << "j = " << j << " multiindex = ("
	     << i0 << ","<< i1 << ","<< i2 << ")"
	     << endl;
      }
      
      // cases:
      // a. corners
      if(i0==nSub) {
	facePoints[f].push_back( face->halfEdges[0]->vertex->id );
	if(debug) cout << "Corner 0: " << facePoints[f].back() << endl;
      } else if(i1==nSub) {
	facePoints[f].push_back( face->halfEdges[1]->vertex->id );
	if(debug) cout << "Corner 1: " << facePoints[f].back() << endl;
      } else if(i2==nSub) {
	facePoints[f].push_back( face->halfEdges[2]->vertex->id );
	if(debug) cout << "Corner 2: " << facePoints[f].back() << endl;
      }
      // b. edges
      else if(i0==0) { 
	int id = face->halfEdges[2]->id;
	facePoints[f].push_back( halfEdgePoints[id][i2-1] );
	if(debug) cout << "Edge i0=0: " << facePoints[f].back() << endl;
      } else if(i1==0) { 
	int id = face->halfEdges[0]->id;
	facePoints[f].push_back( halfEdgePoints[id][i0-1] );
	if(debug) cout << "Edge i1=0: " << facePoints[f].back() << endl;
      }	else if(i2==0) { 
	int id = face->halfEdges[1]->id;
	facePoints[f].push_back( halfEdgePoints[id][i1-1] );
	if(debug) cout << "Edge i2=0: " << facePoints[f].back() << endl;
      } else {
	// c. interior: linearly interpolate using multiindex
	int v0 = face->halfEdges[0]->vertex->id;
	int v1 = face->halfEdges[1]->vertex->id;
	int v2 = face->halfEdges[2]->vertex->id;
	tvmet::Vector<double, 3> p;
	p = i0*points[v0] + i1*points[v1] + i2*points[v2];
	p /= nSub;
	facePoints[f].push_back( npts++ );
	points.push_back(p);

	if(debug) {
	  cout << "Interior (" << i0 << ","<< i1 << ","<< i2 << "): " 
	       << facePoints[f].back() << endl;
	}
      }
      // increment multiindex
      if(i1>0)      {i1--;i2++;}
      else if(i0>0) {i0--;i1=nSub-i0;i2=0;}
      else          {i2++;}

    }
    if(debug){
      cout << "Added " << facePoints[f].size() << " points to face " << f
	   << endl;
    }

  }

  // 4. build new connectivities
  int ntriNew = ntriOld*nSub*nSub;
  HalfEdgeMesh::ConnectivityContainer connectivitiesNew;//( ntriNew );

  if(debug) {
    cout << endl << "Building connectivities for " << ntriNew << " new faces."
	 << endl;
  }
  for(int f=0; f<nF; f++) {
    if(debug) cout << endl << "Face " << f << endl;
    int i0=nSub, i1=0, i2=0;
    while( i2 < nSub && i1 < nSub ) {
      // j2 - j3
      // |  \  |
      // j0 - j1
      int j0 = i2+((i1+i2)*(i1+i2+1))/2;
      int j1 = i2+((i1+1+i2)*(i1+i2+2))/2;
      int j2 = i2+1+((i1+i2+1)*(i1+i2+2))/2;
      int j3 = i2+1+((i1+i2+2)*(i1+i2+3))/2;
      if(debug) {
	cout << "multiindex = ("
	     << i0 << ","<< i1 << ","<< i2 << ")"
	     << endl
	     << " j0 = " << j0
	     << " j1 = " << j1
	     << " j2 = " << j2
	     << " j3 = " << j3
	     << endl;
      }
      HalfEdgeMesh::TriangleConnectivity c;
      c = facePoints[f][j0], facePoints[f][j1], facePoints[f][j2];
      connectivitiesNew.push_back(c);
      if(debug) {
	cout << "New face " << connectivitiesNew.size() 
	     <<" (" << c(0) << ","<< c(1) << ","<< c(2) 
	     << ") from old face " << f 
	     <<" (" << connectivitiesOld[f](0) << ","
	     << connectivitiesOld[f](1) << ","
	     << connectivitiesOld[f](2) << ")"
	     << endl;
      }
      if(j3 < ptsPerFace) {
	c = facePoints[f][j1], facePoints[f][j3], facePoints[f][j2];
	connectivitiesNew.push_back(c);
	if(debug) {
	  cout << "New face " << connectivitiesNew.size() 
	       <<" (" << c(0) << ","<< c(1) << ","<< c(2) 
	       << ") from old face " << f 
	       <<" (" << connectivitiesOld[f](0) << ","
	       << connectivitiesOld[f](1) << ","
	       << connectivitiesOld[f](2) << ")"
	       << endl;
	}
      }
      // increment multiindex
      if(i1>0)      {i1--;i2++;}
      else if(i0>0) {i0--;i1=nSub-i0;i2=0;}
      else          {i2++;}
    }
  }
  
  cout << "Refined mesh has "
       << npts << " vertices and " << connectivitiesNew.size() << " faces."
       << endl;
  if(debug) {
    cout << "Number of new vertices should be " 
	 << vtkMeshOld->GetNumberOfPoints() + 
      heMeshOld.halfEdges.size()/2*(nSub-1) + 
      nF*( ((nSub+1)*(nSub+2))/2 - 3*nSub )
	 << endl
	 << "Number of new faces should be " <<nSub*nSub*nF 
	 << endl;
  }
    
  // print out new mesh to a vtk file
  char  fileName[50];
  sprintf(fileName, "%s-sub%d.vtk",modelName.c_str(),nSub);
  cout << "Saving refined mesh as " << fileName << endl;

  std::ofstream ofs(fileName);
  if (!ofs) {
    std::cout << "can not open output ("
	      << fileName
	      << ") file." << std::endl;
    exit(0);
  }
  
  // print vertices
    ofs << "# vtk DataFile Version 2.0\n"
	<< modelName << " with edges subdivisded " << nSub << " times."
	<< std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << npts << "  double" << std::endl;
    ofs.setf(std::ios_base::scientific,std::ios_base::floatfield);
    for (int i=0; i<npts; i++) {
      ofs << std::setprecision(16) 
	  << std::setw(24) << points[i](0) 
	  << std::setw(24) << points[i](1) 
	  << std::setw(24) << points[i](2) 
	  << std::endl;
    }

    // print faces
    ofs << "POLYGONS  " << connectivitiesNew.size() << "  "
	<< 4*connectivitiesNew.size() << std::endl;
    for (int f=0; f<connectivitiesNew.size(); f++) {
      ofs << 3 << "  "
	  << std::setw(10) << connectivitiesNew[f](0)
	  << std::setw(10) << connectivitiesNew[f](1)
	  << std::setw(10) << connectivitiesNew[f](2)
	  << std::endl;
    }
  return 0;
}

