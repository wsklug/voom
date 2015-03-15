#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"

#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCell.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void insertValenceInVtk(const std::string fileName, vtkSmartPointer<vtkPolyData> mesh);

int main(int argc, char* argv[]){

  string modelName = argv[1];
  string inputFileName = modelName + ".vtk";
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
  vtkSmartPointer<vtkPolyData> mesh = normals->GetOutput();
  mesh->Update();
  if(ifstream(inputFileName.c_str())){
    insertValenceInVtk(inputFileName,mesh);    
  }
}

///////////////////////// INSERTVALENCEINVTK BEGINS ///////////////////////////
/////////////////////////                           //////////////////////////

//The method insertValenceInVtk() inserts valence information in a vtk file
//The calling method must ensure that 'fileName' exists and is a valid vtk file.
//MUST use the extension '.vtk' in 'fileName'.
//'mesh' is a pointer to a vtkPolyData

void insertValenceInVtk(const std::string fileName, vtkSmartPointer<vtkPolyData> mesh){
  ofstream * appendTo;

  //Check that the file exists
  assert( ifstream(fileName.c_str()) );

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName(fileName.c_str());
  
  vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();  
  ds->Update();
  
  //The following vtkUnsignedIntArray will be used to store the
  //number of CELLS in the mesh that share the POINT denoted by the
  //index of the vector
  vtkSmartPointer<vtkUnsignedIntArray> countPointCells = 
    vtkSmartPointer<vtkUnsignedIntArray>::New();
  countPointCells->SetNumberOfValues(mesh->GetNumberOfPoints());
  countPointCells->SetName("Valence");
  
  //cellIds will be used to temporarily hold the CELLS that use a
  //point specified by a point id.
  vtkSmartPointer<vtkIdList> cellIds = 
    vtkSmartPointer<vtkIdList>::New();
  
  //We will use a vtkPolyDataWriter to write our modified output
  //files that will have Capsomer information as well
  vtkSmartPointer<vtkDataWriter> writer = 
    vtkSmartPointer<vtkDataWriter>::New();
  
  if(ds->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
    vtkSmartPointer<vtkUnstructuredGrid> usg 
      = vtkUnstructuredGrid::SafeDownCast(ds);
    usg->BuildLinks();
    for(int p=0; p<ds->GetNumberOfPoints(); p++){
      usg->GetPointCells(p,cellIds);
      countPointCells->SetValue(p,cellIds->GetNumberOfIds());
      cellIds->Reset();
    }      
  }                  
  else if(ds->GetDataObjectType() == VTK_POLY_DATA){
    vtkSmartPointer<vtkPolyData> pd = vtkPolyData::SafeDownCast(ds);
    pd->BuildLinks();
    for(int p=0; p<ds->GetNumberOfPoints(); p++){
      pd->GetPointCells(p,cellIds);
      countPointCells->SetValue(p,cellIds->GetNumberOfIds());
      cellIds->Reset();
    }	
  }
  ds->GetFieldData()->AddArray(countPointCells);
  appendTo = new ofstream(fileName.c_str(),ofstream::app);
  writer->WriteFieldData(appendTo,ds->GetFieldData());
  appendTo->close();
}
