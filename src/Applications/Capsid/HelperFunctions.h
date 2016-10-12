#ifndef _HELPERFUNCTIONS_H
#define _HELPERFUNCTIONS_H


#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unistd.h>
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

#include "Node.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

namespace voom
{

// Function declarations
void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
			double avgEdgeLen, double percentStrain);
void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
			double avgEdgeLen, std::vector<double> percentStrain);
void insertValenceInVtk(std::vector<std::string> fileNames);

std::vector<double> calcEdgeLenAndStdDev(std::vector<DeformationNode<3>*> a,std::vector< tvmet::Vector<int,3> > b);
 
}
#endif
