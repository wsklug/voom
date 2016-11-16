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

double* getBinLimits(vtkSmartPointer<vtkPolyData> poly,
	vtkSmartPointer<vtkIdList> points);
std::complex<double> getIntegrandVal(int l, int m, double *, double, double);

int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << "Usage: sphHarmCoeff <inFile> <outFile> <l_max>"
			<< endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	string inputFile = argv[1];
	string outputFile = argv[2];
	int l_max = std::atoi(argv[3]);

	//Read the input file
	assert(ifstream(inputFile.c_str()));
	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(inputFile.c_str());
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

	vtkSmartPointer<vtkCellArray> bins = poly->GetPolys();
	vtkSmartPointer<vtkIdList> points
		= vtkSmartPointer<vtkIdList>::New();

	ofstream output;
	output.open(outputFile.c_str());
	output << "#l\t\tm\t\tClm" << std::endl;

	//Set the quad-order and quad-points here
	int quadLen = 3;
	double xi[] = { -0.7746, 0.0, 0.7746 };
	double w[] = { 0.555556, 0.888889, 0.555556 };

	for (int l = 0; l < l_max; l++) {
		std::complex<double> Clm = 0;

		for (int m = -l; m <= l; m++) {

			bins->InitTraversal();
			int index = 0;
			while (bins->GetNextCell(points)) {
				double * binLimits = getBinLimits(poly, points);
				double rho = density->GetTuple1(index++);
				double phi1, phi2, theta1, theta2;
				theta1 = binLimits[0];
				theta2 = binLimits[1];
				phi1 = binLimits[2];
				phi2 = binLimits[3];

				double Jacobian = (phi2 - phi1)*(theta2 - theta1) / 4;

				for (int i = 0; i < quadLen; i++) {
					for (int j = 0; j < quadLen; j++) {
						std::complex<double> currVal =
							getIntegrandVal(l, m, binLimits, xi[i], xi[j]);
						Clm += rho*w[i] * w[j] * currVal*Jacobian;
					}
				}
			}
			output << l << "\t\t" << m << "\t\t" << Clm
				<< std::endl;
		}
	}
	t2 = clock();
	double diff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
	output.close();
	return 0;
}

//Function to get theta and phi limits of a cell
double * getBinLimits(vtkSmartPointer<vtkPolyData> poly,
	vtkSmartPointer<vtkIdList> points) {
	double theta_max = 0, theta_min = M_PI,
		phi_max = 0, phi_min = 2*M_PI;

	for (int i = 0; i < points->GetNumberOfIds(); i++) {

		//Get Cartesian coordinates for each point
		double *xyz = poly->GetPoint(points->GetId(i));

		//Skip the "poles" of the sphere
		if ((std::abs(xyz[0]) < 1e-8) &&
			(std::abs(xyz[1]) < 1e-8)) {
			if (std::abs(xyz[2] - 1) < 1e-8)
				theta_min = 0;
			if (std::abs(xyz[2] + 1) < 1e-8)
				theta_max = M_PI;
			continue;
		}

		//Convert to spherical coordinates (phi,theta)
		double phi = atan2(xyz[1], xyz[0]);
		double theta = acos(xyz[2]);

		phi = (phi < 0) ? (2*M_PI + phi) : phi;

		//Compare to update max and min values
		theta_max = std::max(theta, theta_max);
		theta_min = std::min(theta, theta_min);
		phi_min = std::min(phi, phi_min);
		phi_max = std::max(phi, phi_max);
	}
	//Checking for the last bin along phi direction
	//Assuming we will never make longitudinal bin width
	//larger than 45 degrees which is already too much 
	if ((phi_max - phi_min) > M_PI_4) {
		phi_min = phi_max;
		phi_max = 2*M_PI;
	}

	double binLimits[] = { theta_min, theta_max, phi_min, phi_max };

	return binLimits;

}

std::complex<double> getIntegrandVal(int l, int m,
	double* binLimits, double xi, double eta) {

	double phi1, phi2, theta1, theta2;
	theta1 = binLimits[0];
	theta2 = binLimits[1];
	phi1 = binLimits[2];
	phi2 = binLimits[3];

	double theta, phi;

	theta = theta1 + (xi + 1)*(theta2 - theta1) / 2;
	phi = phi1 + (eta + 1)*(phi2 - phi1) / 2;

	std::complex<double> currVal;

	//We want comple conjugate of Y_l^m = (-1)^m*Y_l^(-m)
	currVal = pow(-1, m)*boost::math::spherical_harmonic(l, -m, theta, phi)*sin(theta);

	return currVal;
}