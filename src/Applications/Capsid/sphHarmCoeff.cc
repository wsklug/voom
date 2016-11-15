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
#include <tvmet/Vector.h>
#include <unistd.h>
#include <limits>
#include <algorithm>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkIdList.h>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;

void cellToPoint(std::string file,
	std::string outputFile);

int main(int argc, char* argv[]) {

	if (argc < 3) {
		cout << "Usage: cellDataToPointData <inFile> <outFile>"
			<< endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	string inputFile = argv[1];
	string outputFile = argv[2];

	cellToPoint(inputFile, outputFile);
	t2 = clock();
	double diff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

	return 0;
}

void cellToPoint(std::string inFile, std::string outFile) {
	//Read the input file
	assert(ifstream(inFile.c_str()));
	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(inFile.c_str());
	reader->ReadAllScalarsOn();
	reader->ReadAllVectorsOn();
	vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
	reader->Update();

	vtkSmartPointer<vtkDataArray> density
		= poly->GetCellData()->GetScalars("Density");

	bool densityExists = false;
	if (density.GetPointer()) {
		densityExists = true;
	}
	else {
		std::cout << "Input file does not have density data."
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	ofstream output;
	output.open(outFile.c_str());

	vtkSmartPointer<vtkCellArray> bins = poly->GetPolys();
	vtkSmartPointer<vtkIdList> points
		= vtkSmartPointer<vtkIdList>::New();
	bins->InitTraversal();
	int index = 0;
	while (bins->GetNextCell(points)) {
		double x_sum = 0, y_sum = 0, z_sum = 0;
		double * xyz;
		int pointCount = points->GetNumberOfIds();
		for (int i = 0; i < pointCount; i++) {
			//Get Cartesian coordinates for each point
			xyz = poly->GetPoint(points->GetId(i));
			x_sum += xyz[0];
			y_sum += xyz[1];
			z_sum += xyz[2];
		}
		x_sum /= pointCount;
		y_sum /= pointCount;
		z_sum /= pointCount;

		tvmet::Vector<double, 3> temp(x_sum, y_sum, z_sum);
		tvmet::Vector<double, 3> centroid;
		centroid = temp / tvmet::norm2(temp);

		//Convert to spherical coordinates (phi,theta)
		double phi = (180 / M_PI)*atan2(centroid(1), centroid(0));
		double theta = (180 / M_PI)*acos(centroid(2));
		phi = (phi < 0) ? (360 + phi) : phi;

		double currDensity = density->GetTuple1(index++);

		output << "{" << theta << "," << phi << "},"
			<< currDensity << "}," << std::endl;
	}

	output.close();
	return;
}