#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define _USE_MATH_DEFINES
#include <cmath>
#endif
#include <string>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include <limits>
#include <algorithm>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include "VoomMath.h"

#include <Eigen/Dense>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;

void densityPlotter(std::vector<std::string> fileNames,
	std::string outputFile,
	int lat_res, int long_res);

int main(int argc, char* argv[]) {

	if (argc < 3) {
		cout << "Usage: reOrient <inpFile> <outFile>" << endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	string inputFile = argv[1];
	string outFile = argv[2];

	//Read the input file
	assert(ifstream(inputFile.c_str()));
	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->ReadAllScalarsOn();
	reader->ReadAllVectorsOn();
	vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
	reader->Update();

	// Prepare Eigen matrices
	Eigen::Matrix3Xd target(3, 12), currPentamers(3, 12), rotatedPenta(3,12),
		current(3, poly->GetNumberOfPoints()), 
		newPoints(3, poly->GetNumberOfPoints());

	unsigned short nCells;
	vtkSmartPointer<vtkIdList> cellIds =
		vtkSmartPointer<vtkIdList>::New();

	int pentIndex = 0;
	for (vtkIdType idx = 0; idx < poly->GetNumberOfPoints(); idx++) {
		double pts[] = { 0,0,0 };
		poly->GetPoint(idx, pts);
		current.col(idx) << pts[0], pts[1], pts[2];
		poly->GetPointCells(idx, cellIds);
		if (cellIds->GetNumberOfIds() == 5) {
			Eigen::Vector3d tempVec(pts[0], pts[1], pts[2]);
			currPentamers.col(pentIndex++) = tempVec/tempVec.norm() 
				+ idx*0.1*Eigen::Vector3d::UnitX();
		}
	}
	
	//Hard-coded values for 5-fold site along z-axis
	target << 0, 0, 0.894427, 0.276393, -0.723607, -0.723607,
		0.276393, 0.723607, -0.276393, -0.894427, -0.276393, 0.723607,
		0, 0, 0, 0.850651, 0.525731, -0.525731,
		-0.850651, 0.525731, 0.850651, 0, -0.850651, -0.525731,
		1, -1, 0.447214, 0.447214, 0.447214, 0.447214, 
		0.447214, -0.447214, -0.447214, -0.447214, -0.447214, -0.447214;

	for (int col = 0; col < 12; col++) {
		double temp = target(0, col) + col*0.1;
		target(0, col) = temp;
	}

	Eigen::Affine3d A;

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> tempPointSet1 = vtkSmartPointer<vtkPoints>::New();
	tempPointSet1->SetNumberOfPoints(12);
	vtkSmartPointer<vtkPoints> tempPointSet2 = vtkSmartPointer<vtkPoints>::New();
	tempPointSet2->SetNumberOfPoints(12);

	for (int z = 0; z < 12; z++) {
		double pts[] = {target(0,z),target(1,z),target(2,z)};
		tempPointSet1->SetPoint(z,pts);
		pts[0] = currPentamers(0, z);
		pts[1] = currPentamers(1, z);
		pts[2] = currPentamers(2, z);
		tempPointSet2->SetPoint(z, pts);
	}

	vtkSmartPointer<vtkPolyDataWriter> writer
		= vtkSmartPointer<vtkPolyDataWriter>::New();
	pd->SetPoints(tempPointSet1);
	writer->SetFileName("Target.vtk");
	writer->SetInputData(pd);
	writer->Update();
	writer->Write();

	pd->SetPoints(tempPointSet1);
	writer->SetFileName("MyPentamers.vtk");
	pd->SetPoints(tempPointSet2);
	writer->Update();
	writer->Write();

		A = voom::Find3DAffineTransform(currPentamers, target);

		for (int col = 0; col < current.cols(); col++) {
			newPoints.col(col) = A.linear()*current.col(col)
				+ A.translation();
		}

		for (int col = 0; col < currPentamers.cols(); col++) {
			Eigen::Vector3d temp = currPentamers.col(col) - col*0.1*Eigen::Vector3d::UnitX();
			rotatedPenta.col(col) = A.linear()*temp
				+ A.translation();
		}

		std::cout << "Translation: " << A.translation() <<std::endl;

	for (int z = 0; z < 12; z++) {
		double pts[] = { rotatedPenta(0,z),rotatedPenta(1,z),
			rotatedPenta(2,z) };
		tempPointSet1->SetPoint(z, pts);
	}

	pd->SetPoints(tempPointSet1);
	writer->SetFileName("RotatedPentamers.vtk");
	writer->Update();
	writer->Write();

	vtkSmartPointer<vtkPoints> newPointSet =
		vtkSmartPointer<vtkPoints>::New();
	newPointSet->SetNumberOfPoints(poly->GetNumberOfPoints());

	for (int index = 0; index < newPoints.cols(); index++) {
		double currPoint[] = { newPoints(0,index),
			newPoints(1,index), newPoints(2, index) };
		newPointSet->SetPoint(index, currPoint);
	}

	poly->SetPoints(newPointSet);
	writer->SetFileName(outFile.c_str());
	writer->SetInputData(poly);
	writer->Update();
	writer->Write();

	t2 = clock();
	double diff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

	return 0;
}
