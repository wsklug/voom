#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <iomanip>
#include <limits>

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

#define PI 3.14159265

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

using namespace std;
using namespace tvmet;

std::vector<std::string> &split(const std::string &s, char delim, 
				std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

int main(int argc, char* argv[])
{
  //Usage: ./rotatePoints <input> <output> <angle in degrees> \
  //                   [-c "id1,id2,id3"] | [-m "id1,id2"]
  //
  //The angle must be entered in degrees
  //
  //! Description of the command line options 
  //-c "id1,id2,id3" : 
  // Calculate the axis of rotation from 
  // centroid of three points with ids id1, id2 
  // and id3 in the VTK file
  //
  //-m "id1,id2" : 
  // Calculate the axis of rotation from 
  // midpoint of two points
  //
  //-a "x,y,z"
  // Enter the axis vector directly
  //


  std::string fileName(argv[1]);
  std::string outputFile(argv[2]);
  double theta = atof(argv[3]);
  //convert angle to radians
  theta *= PI/180.0;

  //Check that the file exists
  assert(ifstream(fileName.c_str()));
  
  vtkSmartPointer<vtkPolyDataReader> reader = 
    vtkSmartPointer<vtkPolyDataReader>::New();
  reader->SetFileName(fileName.c_str());
  reader->ReadAllScalarsOn();
  reader->ReadAllVectorsOn();
  
  vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();  
  reader->Update();
   
  int numPoints = pd->GetNumberOfPoints();
  double cos_t = cos(theta);
  double sin_t = sin(theta);

  vtkSmartPointer<vtkDataArray> displacements 
    = pd->GetPointData()->GetVectors("displacements");

  tvmet::Vector<double,3> axis;

  for(int z=4 ; z < argc; z++ ){

    if(std::string(argv[z]) == "-c" || std::string(argv[z]) == "-m"){

      std::string pointIdsString(argv[z+1]);  
      std::vector<int> pointIds;
      std::vector<string> elems = split(pointIdsString,',');
  
      for(std::vector<string>::iterator i=elems.begin(); 
	  i != elems.end(); i++){
	pointIds.push_back(atoi((*i).c_str()));
      }
      z++;

      // This works for both centroid and midpoint  
      tvmet::Vector<double,3> axis_x,axis_d,temp1,temp2,sum;
      sum = 0.0,0.0,0.0;
      for(std::vector<int>::iterator i= pointIds.begin();
	  i != pointIds.end(); i++){
	pd->GetPoint(*i, &(axis_x[0]));
	displacements->GetTuple(*i,&(axis_d[0]));
	temp1 = sum;
	temp2 = temp1 + axis_x + axis_d;
	sum = temp2;
      }
    
      for(int i=0; i < 3; i++){
	axis[i] = sum[i]/(pointIds.size());
      }      

    }

 else if(std::string(argv[z]) == "-a"){
      std::string axisString(argv[z+1]);
      std::vector<string> elems = split(axisString,',');

      for(int i=0; i<3; i++){
	axis[i] = atof((elems[i]).c_str());
      }
      z++;
    }
  }
  
  tvmet::Vector<double,3> k;
  k = tvmet::normalize(axis);

  std::cout<< "Axis of rotation = "
	   << k << std::endl;
  std::cout<<"Angle of rotation = "
	   << theta << std::endl;

  vtkSmartPointer<vtkDoubleArray> rotatedPoints 
    = vtkSmartPointer<vtkDoubleArray>::New();
  rotatedPoints->SetNumberOfComponents(3);
  rotatedPoints->SetNumberOfTuples(numPoints);

  vtkSmartPointer<vtkDoubleArray> rotatedDisp 
    = vtkSmartPointer<vtkDoubleArray>::New();
  rotatedDisp->SetNumberOfComponents(3);
  rotatedDisp->SetNumberOfTuples(numPoints);
  rotatedDisp->SetName("displacements");

  for(int z=0; z < numPoints; z++){
    tvmet::Vector<double,3> v;
    tvmet::Vector<double,3> v_rot;
    tvmet::Vector<double,3> d;
    tvmet::Vector<double,3> d_rot;
    pd->GetPoint(z, &(v[0]));
    displacements->GetTuple(z,&(d[0]));

    v_rot = v*cos_t + cross(k,v)*sin_t + dot(k,v)*(1-cos_t)*k;
    d_rot = d*cos_t + cross(k,d)*sin_t + dot(k,d)*(1-cos_t)*k;
    
    rotatedPoints->SetTuple(z,&(v_rot[0]));
    rotatedDisp->SetTuple(z,&(d_rot[0]));
  }
  
  vtkSmartPointer<vtkPoints> rotPoints 
    = vtkSmartPointer<vtkPoints>::New();
  rotPoints->SetData(rotatedPoints);

  pd->SetPoints(rotPoints);
  pd->GetPointData()->AddArray(rotatedDisp);
  
  vtkSmartPointer<vtkPolyDataWriter> writer = 
    vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(outputFile.c_str());
  writer->SetInputData(pd);
  writer->Write();     
 
}
