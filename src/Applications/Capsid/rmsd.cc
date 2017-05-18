#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
#define _USE_MATH_DEFINES
#include <cmath>
#endif
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>
#include <algorithm>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

#if VTK_MAJOR_VERSION < 6
# define SetInputData SetInput
#endif

using namespace std;
using namespace tvmet;

int main(int argc, char* argv[]) {
	if (argc < 9) {
		cout << "Usage: rmsd <fileBaseName>"
			<< "<startNum> <step> <endNum>"
			<< "<outputFile> <x> <y> <z>" << std::endl;
		return(0);
	}

	clock_t t1, t2;
	t1 = clock();

	double boundingBox[3] = { 0.0,0.0,0.0 };

	std::string fileBaseName = argv[1];
	int startNum = std::atoi(argv[2]);
	int step = std::atoi(argv[3]);
	int endNum = std::atoi(argv[4]);
	std::string outputFile = argv[5];
	boundingBox[0] = std::atof(argv[6]);
	boundingBox[1] = std::atof(argv[7]);
	boundingBox[2] = std::atof(argv[8]);

	std::stringstream sstm;

	ofstream dataFile;
	dataFile.open(outputFile.c_str());
	dataFile << "#Step\tRMSD" << std::endl;

	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();

	vtkSmartPointer<vtkPolyData> pd;
	vtkSmartPointer<vtkDataArray> displacements;
	vtkSmartPointer<vtkDataArray> boundaryCrossed;

	int numPoints = 0;
	tvmet::Vector<double, 3> xi, xj, dispi, dispj;
	tvmet::Vector< double, 3 > diff(0, 0, 0);
	tvmet::Vector<double, 3> bcci(0, 0, 0), bccj(0,0,0);
	double Li = 0;//Lindemann parameter for particle i
	double a = 1.0;//Assuming unit lattice spacing

	//Identify nearest neighbour
	sstm << fileBaseName << "-" << startNum << ".vtk";
	string fileName = sstm.str();
	sstm.str("");
	sstm.clear();
	reader->SetFileName(fileName.c_str());
	pd = reader->GetOutput();
	reader->Update();
	numPoints = pd->GetNumberOfPoints();
	std::vector<int> neighbour(numPoints, -1);
	for (int p = 0; p < numPoints; p++) {
		pd->GetPoint(p, &(xi(0)));
		double dist_min = 10e6;
		for (int q = 0; q < numPoints; q++) {
			if (p != q) {
				pd->GetPoint(q, &(xj(0)));
				double dist = tvmet::norm2(xi - xj);
				if (dist < dist_min) {
					dist_min = dist;
					neighbour[p] = q;
				}
			}
		}
	}

	for (int f = startNum; f <= endNum; f = f + step) {
		sstm << fileBaseName << "-" << f << ".vtk";
		string fileName = sstm.str();
		sstm.str("");
		sstm.clear();
		reader->SetFileName(fileName.c_str());
		pd = reader->GetOutput();
		reader->Update();
		displacements =
			pd->GetPointData()->GetVectors("displacements");
		boundaryCrossed =
			pd->GetPointData()->GetArray("boundaryCrossCount");
		Li = 0;

		// Main calculation for Lindemann parameter
		for (int i = 0; i < numPoints; i++) {
			int j = neighbour[i];
			displacements->GetTuple(i, &(dispi(0)));
			displacements->GetTuple(j, &(dispj(0)));
			boundaryCrossed->GetTuple(i,&(bcci(0)));
			boundaryCrossed->GetTuple(j,&(bccj(0)));
			for (int k = 0; k < 3; k++) {
				dispi(k) = dispi(k) + bcci(k) * boundingBox[k];
				dispj(k) = dispj(k) + bccj(k) * boundingBox[k];
			}
			diff = dispi - dispj;
			Li += tvmet::dot(diff,diff);
		}
		Li /= (2*numPoints*a*a);
		dataFile << f - startNum << "\t" << Li << std::endl;
	}

	dataFile.close();
	t2 = clock();
	double timeDiff = ((float)t2 - (float)t1);
	std::cout << "Post-processing execution time: " << timeDiff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

	return 0;
}