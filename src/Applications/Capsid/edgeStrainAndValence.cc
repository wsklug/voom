#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkNew.h>
#include <vtkLine.h>
#include <vtkIdList.h>
#include <vtkUnsignedIntArray.h>
#include <vtkDataArray.h>

#include "HelperFunctions.h"

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  /*
    Usage: ./edgeStrainAndValence <base-name> <num of files> <Rshift>
    e.g.: ./edgeStrainAndValence T7-relaxed 548 1.00144.
    'Rshift' is the Morse Potential parameter that gives the distance
    between particles at equilibrium separation
  */

  if( argc < 5 ) {
    cout << "Usage: edgeStrainAndValence <base-name> <num of files>"
	 <<" <Rshift> <percentStrain>" << endl;
    return(0);
  }

  clock_t t1,t2;
  t1=clock();
  
  string baseFileName = argv[1];
  int numFiles = std::atoi(argv[2]);
  double Rshift = std::atof(argv[3]);
  double percentStrain = std::atof(argv[4]);

  std::stringstream sstm;

  std::vector<std::string> allVTKFiles;
  allVTKFiles.reserve(numFiles);
  for(int fileNum=0 ; fileNum < numFiles; fileNum++){    
    sstm << baseFileName <<"-"<<fileNum <<".vtk";
    std::string tempString = sstm.str();
    allVTKFiles.push_back(tempString);
    sstm.str("");
    sstm.clear();
  }
  
  insertValenceInVtk(allVTKFiles);
  writeEdgeStrainVtk(allVTKFiles,Rshift,percentStrain);
  t2=clock();
  double diff = ((float)t2-(float)t1);
  std::cout<<"Post-processing execution time: "<<diff/CLOCKS_PER_SEC
  	   <<" seconds"<<std::endl;

}


