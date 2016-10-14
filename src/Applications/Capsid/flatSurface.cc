#include "voom.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <limits>
#include "Node.h"
#include "Model.h"
#include "Lbfgsb.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>

#include "SpringPotential.h"
#include "PotentialBody.h"
#include "BrownianKick.h"
#include "ViscousRegularizer.h"
#include "HelperFunctions.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

typedef DeformationNode<3> Node; // nickname for mechanics nodes
typedef std::vector< Node* > NodeContainer;
typedef blitz::Array<double, 1> Vector_t;
typedef blitz::Array<int, 1> IntArray;

int main(int argc, char* argv[]) {
	clock_t t1, t2;
	t1 = clock();
	if (argc != 2) {
		cout << "usage: " << argv[0] << " <filename>\n";
		return -1;
	}

#ifdef _OPENMP
	std::cout << "************* PARALLELIZATION USING OPENMP ****************"
		<< endl << endl;
#endif

	string inputFileName = argv[1];

	//Numerical viscosity input parameter
	int viterMax;
	bool rescale = true;

	double dt = 9.76e-4;
	int nameSuffix = 0;

	vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();

	std::stringstream sstm;

	reader->SetFileName(inputFileName.c_str());

	vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
	reader->Update();

	std::cout << "Mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
		<< std::endl;

	// create vector of nodes

	int dof = 0;
	std::vector< NodeBase* > nodes;
	std::vector< DeformationNode<3>* > defNodes;

	vtkSmartPointer<vtkDataArray> displacements
		= mesh->GetPointData()->GetVectors("displacements");

	bool displacementsExist = false;

	if (displacements.GetPointer()) {
		displacementsExist = true;
	}

	// read in points
	for (int a = 0; a < mesh->GetNumberOfPoints(); a++) {
		int id = a;
		DeformationNode<3>::Point x;
		mesh->GetPoint(a, &(x[0]));
		if (displacementsExist) {
			tvmet::Vector< double, 3 > disp(0, 0, 0);
			tvmet::Vector< double, 3 > tempSum(0, 0, 0);
			displacements->GetTuple(a, &(disp[0]));
			tempSum = x + disp;
			x = tempSum;
		}
		NodeBase::DofIndexMap idx(3);
		for (int j = 0; j < 3; j++) idx[j] = dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
		nodes.push_back(n);
		defNodes.push_back(n);
	}

	cout << "Number of nodes: " << nodes.size() << endl;

	// read in triangle connectivities
	vector< tvmet::Vector<int, 3> > connectivities;
	tvmet::Vector<int, 3> c;
	int ntri = mesh->GetNumberOfCells();
	connectivities.reserve(ntri);
	std::cout << "Number of triangles: " << ntri << endl;

	for (int i = 0; i < ntri; i++) {
		assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
		for (int a = 0; a < 3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
		connectivities.push_back(c);
	}

	// Calculate side lengths average and std dev of the 
	//equilateral triangles
	std::vector<double> lengthStat =
		calcEdgeLenAndStdDev(defNodes, connectivities);
	double EdgeLength = lengthStat[0];
	double stdDevEdgeLen = lengthStat[1];
	std::cout << "Starting configuration :" << endl
		<< "   Average triangle edge length = " << std::setprecision(10)
		<< EdgeLength << endl
		<< "   Standard deviation = " << std::setprecision(10)
		<< stdDevEdgeLen << endl;
	std::cout.precision(6);

	if (rescale) {
		// Rescale size of the lattice by the average equilateral edge length
		for (int i = 0; i < defNodes.size(); i++) {
			DeformationNode<3>::Point x;
			x = defNodes[i]->point();
			x *= 1.0 / EdgeLength;
			defNodes[i]->setPoint(x);
			defNodes[i]->setPosition(x);
		}
		//Recalculate edge lengths and dependent quantities
		lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
		EdgeLength = lengthStat[0];
		std::cout << "After rescaling :" << endl
			<< "   Average triangle edge length = " << std::setprecision(10)
			<< EdgeLength << endl
			<< "   Standard deviation = " << std::setprecision(10)
			<< stdDevEdgeLen << endl;
		std::cout.precision(6);
	}

	double R_e = EdgeLength;

	//******************* READ COOLING SCHEDULE from File *******

	std::ifstream coolFile("cooling.dat");
	assert(coolFile);
	std::vector<vector<double> > coolVec;
	double curr_D, currViterMax, currPrintStep, currSpringConstant;

	std::string headerline;
	std::getline(coolFile, headerline);

	while (coolFile >> curr_D >> currViterMax >> currPrintStep
		>> currSpringConstant) {
		std::vector<double> currLine;
		currLine.push_back(curr_D);
		currLine.push_back(currViterMax);
		currLine.push_back(currPrintStep);
		currLine.push_back(currSpringConstant);
		coolVec.push_back(currLine);
	}
	coolFile.close();

	string fname = inputFileName.substr(0, inputFileName.find("."));
	string iName;
	string rName;
	string actualFile;

	ofstream myfile;
	std::string dataOutputFile = "BrownianRelax.dat";

	myfile.open(dataOutputFile.c_str());
	myfile << "#Step" << " ParaviewStep" << "\t" << "DiffusionCoeff"
		<< "\t" << "SpringEnergy" << "\t" << "BrownEnergy" << "\t"
		<< "ViscousEnergy" << "\t" << "Total Functional" << "\t"
		<< "MeanSquareDisp"
		<< endl;

	//Parameters for the l-BFGS solver
	int m = 5;
	int maxIter = 1e5;
	double factr = 1.0e+1;
	double pgtol = 1e-7;
	int iprint = 100;
	Lbfgsb solver(3 * nodes.size(), m, factr, pgtol, iprint, maxIter);

	double PotentialSearchRF = 1.2*R_e;
	double springConstant = 1.0;
	double diffusionCoeff = 0;;
	double Cd = 0;
	double viscosity = 0;
	double bkEnergy;
	double vrEnergy;
	double springEnergy;
	double energy;
	int printStep, stepCount = 0;
	int paraviewStep = -1;

	//Setting bounds for z-axis degree of freedom for all nodes
	IntArray nbd(3 * defNodes.size());
	Vector_t l(3 * defNodes.size()), u(3 * defNodes.size());
	l = 0;
	u = 0;
	nbd = 0;
	for (int nbd_idx = 0; nbd_idx < 3 * defNodes.size(); nbd_idx++) {
		if (nbd_idx % 3 == 2) {
			nbd(nbd_idx) = 2;
		}
	}
	// Set bounds for the solver;
	solver.setBounds(nbd, l, u);

	//***************** MAIN LOOP ****************************//

	for (int q = 0; q < coolVec.size(); q++) {

		diffusionCoeff = coolVec[q][0];
		Cd = 1.0 / diffusionCoeff;
		viterMax = coolVec[q][1];
		printStep = (int)coolVec[q][2];
		springConstant = coolVec[q][3];
		viscosity = Cd / dt;
		PotentialSearchRF = 1.2*R_e;

		std::cout << "Viscosity Input Parameters:" << std::endl
			<< " Cd = " << Cd << std::endl
			<< "  D = " << diffusionCoeff << std::endl
			<< " dt = " << dt << std::endl;
		SpringPotential SpringMat(springConstant, R_e);
		PotentialBody * SpringBody = new
			PotentialBody(&SpringMat, defNodes, PotentialSearchRF);
		ViscousRegularizer vr(SpringBody->nodes(), viscosity);
		SpringBody->pushBack(&vr);
		BrownianKick bk(defNodes, Cd, diffusionCoeff, dt);
		SpringBody->pushBack(&bk);

		//Create Model
		Model::BodyContainer bdc;
		bdc.push_back(SpringBody);
		Model model(bdc, nodes);

		bool checkConsistency = false;
		if (checkConsistency) {
			std::cout << "Checking consistency......" << std::endl;
			bk.update2DKick();
			SpringBody->checkConsistency(true);
		}

		//***************************  INNER SOLUTION LOOP ***************************//  

		for (int viter = 0; viter < viterMax; viter++) {

			std::cout << std::endl
				<< "VISCOUS ITERATION: " << viter + stepCount
				<< "\t viscosity = " << vr.viscosity()
				<< std::endl
				<< std::endl;

			// Impose the Brownian perturbation
			bk.update2DKick();
			SpringBody->recomputeNeighbors(PotentialSearchRF);

			bool printBeforeSolve = true;
			if (printBeforeSolve) {
				//We will print only after every currPrintStep iterations
				if (viter % printStep == 0) {
					sstm << fname << "-initial-" << nameSuffix;
					rName = sstm.str();
					model.print(rName);
					sstm << "-bd1.vtk";
					actualFile = sstm.str();
					sstm.str("");
					sstm.clear();
					sstm << fname << "-initial-" << nameSuffix << ".vtk";
					rName = sstm.str();
					std::rename(actualFile.c_str(), rName.c_str());
					sstm.str("");
					sstm.clear();
				}
			}
			solver.solve(&model);
			bkEnergy = bk.energy();
			vrEnergy = vr.energy();
			springEnergy = SpringBody->energy() - bkEnergy - vrEnergy;
			energy = solver.function();

			std::cout << "ENERGY:" << std::endl
				<< "viscous energy  = " << vrEnergy << std::endl
				<< "Brownian energy = " << bkEnergy << std::endl
				<< "Spring energy   = " << springEnergy << std::endl
				<< "  total energy  = " << energy << std::endl
				<< std::endl;
			std::cout << "VISCOSITY: " << std::endl
				<< "          velocity = " << vr.velocity() << std::endl
				<< " updated viscosity = " << vr.viscosity() << std::endl
				<< std::endl;


			//*********************************************************//

			std::cout << "Shape relaxed." << std::endl
				<< "Energy = " << energy << std::endl;


			//********** Print relaxed configuration ************// 

			//We will print only after every currPrintStep iterations
			if (viter % printStep == 0) {
				paraviewStep++;
				sstm << fname << "-relaxed-" << nameSuffix;
				rName = sstm.str();
				model.print(rName);
				sstm << "-bd1.vtk";
				actualFile = sstm.str();
				sstm.str("");
				sstm.clear();
				sstm << fname << "-relaxed-" << nameSuffix << ".vtk";
				rName = sstm.str();
				std::rename(actualFile.c_str(), rName.c_str());
				sstm.str("");
				sstm.clear();
			}
			//****************************************************//

			double msd = 0; //Mean Squared Displacement
			for (int i = 0; i < defNodes.size(); i++) {
				Vector3D tempDisp;
				tempDisp = (defNodes[i]->point() - defNodes[i]->position());
				msd += tvmet::dot(tempDisp, tempDisp);
			}

			msd /= (R_e*R_e*defNodes.size());

			int paraviewStepPrint;
			paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;

			myfile << nameSuffix++ << "\t\t" << paraviewStepPrint << "\t\t"
				<< diffusionCoeff << "\t\t"
				<< springEnergy << "\t\t"
				<< bkEnergy << "\t\t"
				<< vrEnergy << "\t\t" << energy << "\t\t" << msd
				<< endl;

			// step forward in "time", relaxing viscous energy & forces 
			vr.step();
		}
		delete SpringBody;
	}
	myfile.close();
	t2 = clock();
	float diff((float)t2 - (float)t1);
	std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
}