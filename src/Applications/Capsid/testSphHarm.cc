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
#include <vtkSphereSource.h>
#include <vtkDoubleArray.h>
#include <vtkPolyDataWriter.h>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;

double* getBinLimits(vtkSmartPointer<vtkPolyData> poly,
	vtkSmartPointer<vtkIdList> points);

int main(int argc, char* argv[]) {

	if (argc < 5) {
		cout << "Usage: testSphHarm <lat_res> <long_res> <inFile> <outFile>"
			<< endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	int lat_res = std::atoi(argv[1]);
	int long_res = std::atoi(argv[2]);
	string inFile = argv[3];
	string outFile = argv[4];

	std::ifstream input(inFile.c_str());
	int l_temp, m_temp;
	double Clm_temp1, Clm_temp2;

	std::vector<int> l, m;
	std::vector<std::complex<double> > Clm ;

	while ( input >> l_temp >> m_temp >> Clm_temp1 >> Clm_temp2) {
		l.push_back(l_temp);
		m.push_back(m_temp);
		std::complex<double> Clm_temp;
		Clm_temp.real(Clm_temp1);
		Clm_temp.imag(Clm_temp2);
		Clm.push_back(Clm_temp);
	}

	vtkSmartPointer<vtkSphereSource> sp =
		vtkSmartPointer<vtkSphereSource>::New();
	sp->SetRadius(1.0);
	sp->SetThetaResolution(lat_res);
	sp->SetPhiResolution(long_res);
	sp->LatLongTessellationOn();

	vtkSmartPointer<vtkPolyData> poly =
		sp->GetOutput();
	sp->Update();

	vtkSmartPointer<vtkCellArray> bins = poly->GetPolys();
	vtkSmartPointer<vtkIdList> points
		= vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkDoubleArray> scalarDensity =
		vtkSmartPointer<vtkDoubleArray>::New();	
	scalarDensity->SetNumberOfComponents(1);
	scalarDensity->SetNumberOfTuples(poly->GetNumberOfCells());
	scalarDensity->SetName("Density");
	scalarDensity->FillComponent(0, 0.0);

	vtkSmartPointer<vtkDoubleArray> vectorDensity =
		vtkSmartPointer<vtkDoubleArray>::New();
	vectorDensity->SetNumberOfComponents(2);
	vectorDensity->SetNumberOfTuples(poly->GetNumberOfCells());
	vectorDensity->SetName("VectorDensity");
	vectorDensity->FillComponent(0, 0.0);
	vectorDensity->FillComponent(1, 0.0);

	bins->InitTraversal();
	int index = 0;
	while (bins->GetNextCell(points)) {
		double * binLimits = getBinLimits(poly, points);
		double phi1, phi2, theta1, theta2;
		theta1 = binLimits[0];
		theta2 = binLimits[1];
		phi1 = binLimits[2];
		phi2 = binLimits[3];

		double theta = (theta1 + theta2) / 2;
		double phi = (phi1 + phi2) / 2;

		std::complex<double> den = 0;

		for (int idx = 0; idx < Clm.size(); idx++) {
			den += Clm[idx]*boost::math::spherical_harmonic(l[idx], 
				m[idx], theta, phi);
		}

		vectorDensity->SetTuple2(index, std::real(den), std::imag(den));
		scalarDensity->SetTuple1(index++,std::real(den));
	}

	poly->GetCellData()->AddArray(scalarDensity);
	poly->GetCellData()->AddArray(vectorDensity);

	vtkSmartPointer<vtkPolyDataWriter> wr =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	wr->SetFileName(outFile.c_str());
	wr->SetInputData(poly);
	wr->Update();
	wr->Write();

	t2 = clock();
	double diff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
	return 0;
}

//Function to get theta and phi limits of a cell
double * getBinLimits(vtkSmartPointer<vtkPolyData> poly,
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

	double binLimits[] = { theta_min, theta_max, phi_min, phi_max };

	return binLimits;

}
