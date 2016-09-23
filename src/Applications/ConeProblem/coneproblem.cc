#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Model.h"

#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetWriter.h>
#include <vtkPolyDataReader.h>

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
 if (argc != 2)
     {
     cout << "Usage: " << argv[0]
          << " InputPolyDataFile" << endl;
     return EXIT_FAILURE;
     }
 
   //Read the file
   vtkSmartPointer<vtkPolyDataReader> reader =
     vtkSmartPointer<vtkPolyDataReader>::New();
   string inputPolyDatafile = argv[1];
   inputPolyDatafile = inputPolyDatafile + ".vtk";
   reader->SetFileName (inputPolyDatafile.c_str());
   reader->Update();

   vtkSmartPointer<vtkPolyData> polyData =
     reader->GetOutput();
   
   double bounds[6];
   polyData->GetBounds(bounds);
   
   std::cout  << "xmin: " << bounds[0] << " " 
	      << "xmax: " << bounds[1] << std::endl
	      << "ymin: " << bounds[2] << " " 
	      << "ymax: " << bounds[3] << std::endl
	      << "zmin: " << bounds[4] << " " 
	      << "zmax: " << bounds[5] << std::endl;
   
   // Now we have to add 100 points from xmin to xmax and add it to
   // the vtk file

   double xmin = bounds[0];
   double xmax = bounds[1];
   
   int numOfPointsAlready = polyData->GetNumberOfPoints();
   int numPointsToInsert = 100;
   double dx = (xmax - xmin)/(numPointsToInsert-1);
   double x;
   double insertAtId = numOfPointsAlready;
   
   vtkSmartPointer<vtkPoints> points = polyData->GetPoints();

   for(int z=1;z<=numPointsToInsert;z++){
     x = xmin+(z-1)*dx;
     points->InsertPoint(insertAtId++,x,0.0,0.0);
   }
   
   polyData->SetPoints(points);

   // Clean the polydata. This will remove duplicate
   // points that may be // present in the input data.
   vtkSmartPointer<vtkCleanPolyData> cleaner =
   vtkSmartPointer<vtkCleanPolyData>::New();
   cleaner->SetInput (polyData);
   cleaner->Update();
   // Generate a tetrahedral mesh from the input points. By
   // default, the generated volume is the convex hull of the points.
   vtkSmartPointer<vtkDelaunay3D> delaunay3D =
     vtkSmartPointer<vtkDelaunay3D>::New();
   delaunay3D->SetInputConnection (cleaner->GetOutputPort());

   // Write the mesh as an unstructured grid
   vtkSmartPointer<vtkDataSetWriter> writer =
     vtkSmartPointer<vtkDataSetWriter>::New();
   string meshUSG =  inputPolyDatafile + "_mesh.vtk";

   writer->SetFileName (meshUSG.c_str());
   writer->SetInputConnection ( delaunay3D->GetOutputPort() );
   writer->Write();

   return EXIT_SUCCESS;
 
}
