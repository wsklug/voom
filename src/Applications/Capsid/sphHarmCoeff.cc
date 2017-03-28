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

#include <limits>
#include <algorithm>
#include <boost/math/special_functions/spherical_harmonic.hpp>

#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;

std::vector<double> getBinLimits(vtkSmartPointer<vtkPolyData> poly,
	vtkSmartPointer<vtkIdList> points);

 void renormalizeDensity(vtkSmartPointer<vtkPolyData> poly);

int main(int argc, char* argv[]) {

	if (argc < 6) {
		cout << "Usage: sphHarmCoeff <inFile> <outFile> <l_max> <vectorDensityFlag> <renormalizeDensityFlag>"
			<< endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	string inputFile = argv[1];
	string outputFile = argv[2];
	int l_max = std::atoi(argv[3]);
	bool vectorDensityFlag = std::atoi(argv[4]);
	bool renormalizeDensityFlag = std::atoi(argv[5]);

	//Read the input file
	assert(ifstream(inputFile.c_str()));
	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->ReadAllScalarsOn();
	reader->ReadAllVectorsOn();
	vtkSmartPointer<vtkPolyData> poly = reader->GetOutput();
	reader->Update();

	vtkSmartPointer<vtkDataArray> density;

	if (vectorDensityFlag) {
		density = poly->GetCellData()->GetVectors("VectorDensity");
	}
	else {
		density = poly->GetCellData()->GetScalars("Density");
	}

	bool densityExists = false;
	if (density.GetPointer()) {
		densityExists = true;
	}
	else {
		std::cout << "Input file does not have density data."
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	if (renormalizeDensityFlag) {
		renormalizeDensity(poly);
	}

	vtkSmartPointer<vtkCellArray> bins = poly->GetPolys();
	vtkSmartPointer<vtkIdList> points
		= vtkSmartPointer<vtkIdList>::New();

	ofstream output;
	output.open(outputFile.c_str());
	output << "#l\t\tm\t\tClm_real\t\tClm_imag" << std::endl;

	//Set the quad-order and quad-points here
	//3-point quadrature
	/*int quadLen = 3;
	double xi[] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
	double w[] = { 0.5555555555555556,
		0.888888888888889, 0.5555555555555556 };*/

	//4-point quadrature
	int quadLen = 4;
	double quad1 = 0.339981043584856;
	double quad2 = 0.861136311594053;
	double w1 = 0.652145154862546;
	double w2 = 0.347854845137454;

	double xi[] = { -quad2,-quad1 ,quad1 ,quad2 };
	double w[] = { w2,w1,w1,w2 };

	for (int l = 0; l < l_max; l++) {
		for (int m = -l; m <= l; m++) {

			std::complex<double> Clm = 0;
			bins->InitTraversal();
			int index = 0;
			while (bins->GetNextCell(points)) {
			std::vector<double> binLimits = getBinLimits(poly, points);
				std::complex<double> rho;
				if (vectorDensityFlag) {
					double* rhoTemp = density->GetTuple2(index++);
					rho.real(rhoTemp[0]);
					rho.imag(rhoTemp[1]);
				}
				else {
					rho.real(density->GetTuple1(index++));
					rho.imag(0);
				}

				double phi1, phi2, theta1, theta2;
				theta1 = binLimits[0];
				theta2 = binLimits[1];
				phi1 = binLimits[2];
				phi2 = binLimits[3];

				double Jacobian = (phi2 - phi1)*(theta2 - theta1) / 4;

				for (int i = 0; i < quadLen; i++) {
					for (int j = 0; j < quadLen; j++) {
						double theta = theta1 + (xi[i] + 1)*(theta2 - theta1) / 2;
						double phi = phi1 + (xi[j] + 1)*(phi2 - phi1) / 2;
						std::complex<double> currVal =
							pow(-1, m)*boost::math::spherical_harmonic(l, -m, theta, phi)
							*sin(theta);
						Clm += rho * w[i] * w[j] * currVal * Jacobian;
						//Clm += rho * w[i] * w[j] * Jacobian * sin(theta);
					}
				}
			}
			output << l << "\t\t" << m << "\t\t" << std::real(Clm)
				<< "\t\t" << std::imag(Clm) << "\t\t" << std::endl;
		}
	}

	t2 = clock();
	double diff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
	//output.close();
	return 0;
}

//Function to get theta and phi limits of a cell
std::vector<double> getBinLimits(vtkSmartPointer<vtkPolyData> poly,
	vtkSmartPointer<vtkIdList> points) {
	double theta_max = 0, theta_min = M_PI,
		phi_max = 0, phi_min = 2 * M_PI;

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

		phi = (phi < 0) ? (2 * M_PI + phi) : phi;

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
		phi_max = 2 * M_PI;
	}

	std::vector<double> binLimits 
		= { theta_min, theta_max, phi_min, phi_max };

	return binLimits;

}

void renormalizeDensity(vtkSmartPointer<vtkPolyData> poly) {
	vtkSmartPointer<vtkDataArray> vecDensity, density;
		vecDensity = poly->GetCellData()->GetVectors("VectorDensity");
		density = poly->GetCellData()->GetScalars("Density");

	bool densityExists = false;
	bool vecDensityExists = false;
	if (density.GetPointer()) {
		densityExists = true;
	}
	if (vecDensity.GetPointer()) {
		vecDensityExists = true;
	}
	if(!(vecDensityExists | densityExists)) {
		std::cout << "Input file does not have any density data."
			<< std::endl;
		exit(EXIT_FAILURE);
	}
	
	double densityAvg = 0.0;
	std::complex<double> vecDensityAvg = 0.0;

	for (int idx = 0; idx < poly->GetNumberOfCells(); idx++) {
		if (densityExists) {
			densityAvg += density->GetTuple1(idx);
		}
		if (vecDensityExists) {
			double* temp = vecDensity->GetTuple2(idx);
			std::complex<double> tempDen;
			tempDen.real(temp[0]);
			tempDen.imag(temp[1]);
			vecDensityAvg += tempDen;
		}
	}
	densityAvg /= poly->GetNumberOfCells();
	vecDensityAvg /= poly->GetNumberOfCells();

	vtkSmartPointer<vtkDoubleArray> newDensity =
		vtkSmartPointer<vtkDoubleArray>::New();
	newDensity->SetNumberOfComponents(1);
	newDensity->SetNumberOfTuples(poly->GetNumberOfCells());
	newDensity->SetName("Density");
	newDensity->FillComponent(0, 0.0);

	vtkSmartPointer<vtkDoubleArray> newVecDensity =
		vtkSmartPointer<vtkDoubleArray>::New();
	newVecDensity->SetNumberOfComponents(2);
	newVecDensity->SetNumberOfTuples(poly->GetNumberOfCells());
	newVecDensity->SetName("VectorDensity");
	newVecDensity->FillComponent(0, 0.0);
	newVecDensity->FillComponent(1, 0.0);

	for (int idx = 0; idx < poly->GetNumberOfCells(); idx++) {
		if (densityExists) {
			double tempDensity 
				= density->GetTuple1(idx) - densityAvg;
			newDensity->SetTuple1(idx,tempDensity);
		}
		if (vecDensityExists) {
			double* temp = vecDensity->GetTuple2(idx);
			double temp1, temp2;
			temp1 = temp[0] - vecDensityAvg.real();
			temp2 = temp[1] - vecDensityAvg.imag();
			newVecDensity->SetTuple2(idx, temp1, temp2);
		}
	}
	poly->GetCellData()->AddArray(newDensity);
	poly->GetCellData()->AddArray(newVecDensity);

	vtkSmartPointer<vtkPolyDataWriter> wr =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	wr->SetFileName("Renormalized.vtk");
	wr->SetInputData(poly);
	wr->Update();
	wr->Write();

	return;
}
