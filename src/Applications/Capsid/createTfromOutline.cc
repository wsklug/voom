#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <iomanip>
#include <limits>
#include <algorithm>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkLine.h>
#include <vtkIdList.h>
#include <vtkUnsignedIntArray.h>
#include <vtkDataArray.h>
#include <vtkDelaunay3D.h>
#include <vtkCleanPolyData.h>
#include <vtkDataSetSurfaceFilter.h>

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

using namespace std;
using namespace tvmet;

/*
  Usage: ./createTfromOutline <input vtk file name> <output vtk file name>
*/

int main(int argc, char* argv[])
{
  if(argc < 3){
    
    std::cerr << "Usage: "<< argv[0] 
	      << "<input VTK file> <output VTK filename>" << std::endl;
    return 1;
    
  }

  string fileName(argv[1]);
  string outputFile(argv[2]);
  
  //Check that the input file exists
  assert(ifstream(fileName.c_str()));

  vtkSmartPointer<vtkPolyDataReader> reader = 
    vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fileName.c_str());
  
  vtkSmartPointer<vtkPolyDataNormals> normals = 
    vtkSmartPointer<vtkPolyDataNormals>::New();
  normals->SetInputConnection(reader->GetOutputPort());
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  
  vtkSmartPointer<vtkCleanPolyData> cleaner
    = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputConnection(normals->GetOutputPort());

  vtkSmartPointer<vtkPolyData> pd = cleaner->GetOutput();
  cleaner->Update();
 
  int numPoints = pd->GetNumberOfPoints();
  
  std::vector<int> fiveFoldPoints;

  pd->BuildLinks();
  pd->BuildCells();

  //*******************************************************************//
  //Keep only those points that are shared by 2 or 3 hexagons. Discard
  //other points
  //*******************************************************************//
  for(int p=0; p < numPoints; p++){

    vtkSmartPointer<vtkIdList> cellIds = 
      vtkSmartPointer<vtkIdList>::New();
      
    pd->GetPointCells(p,cellIds);
      
    if ( cellIds->GetNumberOfIds() == 2){
      fiveFoldPoints.push_back(p);
    }
  }
  //------------------------------------------------------------------//

  //******************************************************************//
  //Identify the pentagon connectivity
  //******************************************************************//
  int pLSize = fiveFoldPoints.size()/5;
  std::cout << "Number of pentagons = " << pLSize << std::endl;

  std::vector<vector<int> > pentagonList;

  while(fiveFoldPoints.size() > 0){
    
    int currPoint;    
    std::vector<int> pentagon;
    
    currPoint = fiveFoldPoints.back();
    pentagon.push_back(currPoint);
    int index = 0;

    while(pentagon.size() < 5){

      vtkSmartPointer<vtkIdList> cellIds = 
	vtkSmartPointer<vtkIdList>::New();
      
      pd->GetPointCells(pentagon[index++], cellIds);
	
      for(int q=0; q < 2; q++){

	vtkSmartPointer<vtkIdList> cellPoints = 
	  vtkSmartPointer<vtkIdList>::New();
	pd->GetCellPoints( cellIds->GetId(q), cellPoints );

	std::vector<int>::iterator it1;
	std::vector<int>::iterator it2;

	for(int r=0; r < 6; r++){
	  
	  int cellPoint = cellPoints->GetId(r);
	  it1 = find ( fiveFoldPoints.begin(), fiveFoldPoints.end(), cellPoint);
	 
	  if( it1!=fiveFoldPoints.end() ){
	    it2 = find ( pentagon.begin(), pentagon.end(), cellPoint);	  
	    if(it2 == pentagon.end()){	    
	      pentagon.push_back(cellPoint);
	    }    
	  }
	}
      }
    }
    //Need to re-arrange the point ids as per connectivity
    int tempId = pentagon[2];
    pentagon[2] = pentagon[3];
    pentagon[3] = tempId;
    tempId = pentagon[4];
    pentagon[4] = pentagon[3];
    pentagon[3] = tempId;
    
    pentagonList.push_back(pentagon);
    
    std::vector<int>::iterator it3;
    std::vector<int>::iterator it4;
    for( it3=pentagon.begin(); it3!=pentagon.end(); ++it3 ){
      it4 = find ( fiveFoldPoints.begin(), fiveFoldPoints.end(), *it3);
      fiveFoldPoints.erase(it4);
    }
  }
  //-----------------------------------------------------------------------//

  //**************************************************************************//
  // Calculate and store centroid of all hexagons and pentagons.             //
  //*************************************************************************//
  
  std::cout<< "Pentagon points identified" << std::endl;
  
  vtkCellArray* hexagons = pd->GetPolys();

  int numHex = hexagons->GetNumberOfCells();
  
  vtkDoubleArray* finalPoints = vtkDoubleArray::New();
  finalPoints->SetNumberOfComponents(3);
  finalPoints->SetNumberOfTuples( numHex + pLSize);
  int finalIndex = 0;

  vtkIdType npts;
  vtkIdType* pts;
  for(hexagons->InitTraversal();hexagons->GetNextCell( npts, pts);){

    double newPoint[3] = {0,0,0};
    double hexPoint[3];

    for(int s=0; s < npts; s++){

      pd->GetPoint(pts[s], hexPoint);
      newPoint[0] += hexPoint[0];
      newPoint[1] += hexPoint[1];
      newPoint[2] += hexPoint[2]; 
     
    }

    newPoint[0] = newPoint[0]/npts;
    newPoint[1] = newPoint[1]/npts;
    newPoint[2] = newPoint[2]/npts;

    finalPoints->SetTuple(finalIndex++, newPoint );
  }

  for(std::vector< vector<int> >::iterator z = pentagonList.begin(); 
      z != pentagonList.end(); ++z){

    double newPoint[3] = {0,0,0};
    double pentaPoint[3];

    for(std::vector<int>::iterator m = (*z).begin();
	m != (*z).end(); ++m){
     
      pd->GetPoint((*m), pentaPoint);
      newPoint[0] += pentaPoint[0];
      newPoint[1] += pentaPoint[1];
      newPoint[2] += pentaPoint[2];
   
    }

    newPoint[0] = newPoint[0]/5;
    newPoint[1] = newPoint[1]/5;
    newPoint[2] = newPoint[2]/5;

    finalPoints->SetTuple(finalIndex++, newPoint );

  }
  
  std::cout<< "Total number of points = "<< finalIndex << std::endl;  

  //------------------------------------------------------------------------//

  //***********************************************************************//
  //                  Write a new VTK file with finalPoints                //
  //***********************************************************************//
  
  vtkSmartPointer<vtkPoints> newTPoints
    = vtkSmartPointer<vtkPoints>::New();
  newTPoints->SetData(finalPoints);

  vtkSmartPointer<vtkPolyData> newPolyData 
    = vtkSmartPointer<vtkPolyData>::New();
  newPolyData->SetPoints(newTPoints);
  
  vtkSmartPointer<vtkPolyDataWriter> writer = 
    vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(outputFile.c_str());
  writer->SetInputData(newPolyData);
  writer->Update();
  writer->Write();

  //-----------------------------------------------------------------------//
}
