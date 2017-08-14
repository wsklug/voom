/*
 * OPSAsphericity.cc
 *
 *  Created on: Aug 7, 2017
 *      Author: amit
 */
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>
#include "Node.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "OPSBody.h"
#include "BrownianKickOPS.h"
#include "ViscousRegularizer.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include "HelperFunctions.h"

#include <Eigen/Dense>

using namespace tvmet;
using namespace std;
using namespace voom;

typedef tvmet::Vector<double,6> Vector6D;
typedef Eigen::Matrix<double,6,1> Vector6d;
typedef Eigen::Matrix<double,6,Eigen::Dynamic> Matrix6Xd;
typedef Eigen::Transform<double,6,Eigen::Affine> Affine6d;

int main(int argc, char* argv[])
{
	clock_t t1, t2, t3;
	t1 = clock();
	if (argc != 2) {
		cout << "usage: " << argv[0] << " <filename>\n";
		return -1;
	}

	string inputFileName = argv[1];
	string fname = inputFileName.substr(0, inputFileName.find("."));
	string rName;
	std::stringstream sstm;

	//For Morse material
	double epsilon, re, s, K, a, b, aM, aP, aN, aC;
	double percentStrain = 10;
	bool projectOnSphere = false;
	double initialSearchRad = 1.0, finalSearchRad = 1.2;
	int lat_res, long_res, viterMax = 1000;
	int nameSuffix = 0;
	double dt = 9.76e-4;

	//Read epsilon and percentStrain from input file. percentStrain is
	//calculated so as to set the inflection point of Morse potential
	//at a fixed distance relative to the equilibrium separation
	//e.g. 1.1*R_eq, 1.5*R_eq etc.
	std::ifstream miscInpFile("miscInp.dat");
	assert(miscInpFile);
	string temp;
	miscInpFile
	>> temp >> epsilon
	>> temp >> re
	>> temp >> percentStrain
	>> temp >> K
	>> temp >> a
	>> temp >> b
	>> temp >> initialSearchRad
	>> temp >> finalSearchRad
	>> temp >> lat_res
	>> temp >> long_res;
	miscInpFile.close();

	s = (100 / (re*percentStrain))*log(2.0);
	struct OPSParams props = {aM, aP, aN, aC, epsilon, re, s, K, a, b};

	vtkSmartPointer<vtkPolyDataReader> reader =
			vtkSmartPointer<vtkPolyDataReader>::New();
	vtkSmartPointer<vtkPolyData> mesh;

	reader->SetFileName(inputFileName.c_str());
	reader->Update();
	mesh = reader->GetOutput();
	std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
  	    												<< std::endl;

	// create vector of nodes
	int dof = 0;
	std::vector<OPSNode* > nodes;
	std::vector<NodeBase*> baseNodes;
	double Ravg = 0.0;

	// read in points
	for (int i = 0; i < mesh->GetNumberOfPoints(); i++) {
		int id = i;
		NodeBase::DofIndexMap idx(6);
		Vector3D X;
		OPSNode* n;
		for (int j = 0; j < 6; j++) idx[j] = dof++;
		mesh->GetPoint(i, &(X[0]));
		n = new OPSNode(id, idx, X);
		Ravg += tvmet::norm2(X);
		nodes.push_back(n);
		baseNodes.push_back(n);
	}
	assert(nodes.size() != 0);
	Ravg /= nodes.size();
	cout << "Number of nodes: " << nodes.size() << endl
			<< "Initial radius: " << Ravg << endl;

	OPSBody* bd = new OPSBody( nodes, props, initialSearchRad );

	// Calculate side lengths average of the imaginary equilateral triangles
	double EdgeLength = bd->getAverageEdgeLength();

	// Rescale size of the capsid by the average equilateral edge length
	for (int i = 0; i < nodes.size(); i++) {
		Vector3D X;
		X = nodes[i]->referencePosition();
		X *= 1.0 / EdgeLength;
		nodes[i]->setReferencePosAndRotVec(X);
		nodes[i]->setDeformedPosAndRotVec(X);
	}
	//Recalculate edge lengths and capsid radius
	bd->updateSearchRadius(finalSearchRad);
	EdgeLength = bd->getAverageEdgeLength();
	bd->updatePolyDataAndKdTree();
	bd->updateNeighbors();

	Ravg = 0.0;
	for (int i = 0; i < nodes.size(); i++) {
		Vector3D x;
		x = nodes[i]->deformedPosition();
		double tempRadius = tvmet::norm2(x);
		Ravg += tempRadius;
	}
	Ravg /= nodes.size();

	std::cout << "Radius of capsid after rescaling = " << Ravg << endl;

	// Prepare Eigen matrices
	Matrix6Xd initial(6, nodes.size()), current(6, nodes.size());

	//Fill in matrix of the initial state for Kabsch algorithm
	for (int col = 0; col < nodes.size(); col++) {
		Vector3D coord = nodes[col]->deformedPosition();
		for (int row = 0; row < 6; row++) {
			initial(row, col) = coord(row);
		}
	}

	//Prepare output data file
	std::string dataOutputFile;
	sstm << fname << "-output.dat";
	dataOutputFile = sstm.str();
	sstm.str("");
	sstm.clear();

	ofstream myfile;
	myfile.open(dataOutputFile.c_str());
	myfile << "#Step" << "\t"
			<<"ParaviewStep" << "\t"
			<< "DiffusionCoeff" << "\t"
			<< "alphaM" << "\t"
			<< "Asphericity" << "\t"
			//<< "LoopAsphericity" << "\t"
			<< "Radius" << "\t"
			<< "MorseEnergy" << "\t"
			<< "PlanarityEn" << "\t"
			<< "NormalityEn" << "\t"
			<< "CircularityEn" << "\t"
			<< "BrownianEnergy" << "\t"
			<< "ViscosityEnergy" << "\t"
			<< "TotalFunctional"
			<< std::endl;

	// Update the Morse parameters
	s = (100 / (EdgeLength*percentStrain))*log(2.0);
	bd->updateProperty( OPSBody::r, EdgeLength );
	bd->updateProperty( OPSBody::sv, s );

	//Create a l-BFGS-b solver
	int m = 7;
	int maxIter = 1e5;
	double factr = 10.0;
	double pgtol = 1e-7;
	int iprint = 2000;
	Lbfgsb solver(6 * nodes.size(), m, factr, pgtol, iprint, maxIter, true);

	double diffusionCoeff = 4.0*EdgeLength*EdgeLength;
	double Cd = 1.0/diffusionCoeff;
	double viscosity = Cd/dt;

	// Create BrownianKick element
	BrownianKickOPS *bk = new BrownianKickOPS(nodes,Cd,diffusionCoeff,dt);
	bd->addElement( bk );

	// Create ViscousRegularizer element
	ViscousRegularizer *vr = new ViscousRegularizer(baseNodes,viscosity);
	bd->addElement( vr );

	//Create Model
	Model::BodyContainer bdc;
	bdc.push_back(bd);
	Model model(bdc, baseNodes);

	int printStep, stepCount = 0;
	int paraviewStep = -1;

	vtkSmartPointer<vtkSphereSource> sp =
			vtkSmartPointer<vtkSphereSource>::New();
	sp->SetRadius(1.0);
	sp->SetThetaResolution(lat_res);
	sp->SetPhiResolution(long_res);
	sp->LatLongTessellationOn();

	vtkSmartPointer<vtkPolyData> pd =
			sp->GetOutput();
	sp->Update();

	/*
	 * In the next few lines we will identify spherical co-ordinate
	 * limits for each cell in the sphere
	 */
	std::vector<vector<double> > cellLimits = getSphCellLimits(pd, long_res);

	vtkSmartPointer<vtkDoubleArray> binDensity =
			vtkSmartPointer<vtkDoubleArray>::New();
	binDensity->SetNumberOfComponents(1);
	binDensity->SetNumberOfTuples(pd->GetNumberOfCells());
	binDensity->SetName("Density");

	//******************* READ COOLING SCHEDULE from File *******

	std::ifstream coolFile("cooling.dat");
	assert(coolFile);
	std::vector<vector<double> > coolVec;
	double curr_D, currAm, currPercentStrain, currViterMax, currPrintStep;

	std::string headerline;
	std::getline(coolFile, headerline);

	while (coolFile >> curr_D >> currViterMax >> currAm >>
			currPercentStrain >> currPrintStep) {
		std::vector<double> currLine;
		currLine.push_back(curr_D);
		currLine.push_back(currViterMax);
		currLine.push_back(currAm);
		currLine.push_back(currPercentStrain);
		currLine.push_back(currPrintStep);
		coolVec.push_back(currLine);
	}
	coolFile.close();

	//***************************  SOLUTION LOOP ***************************
	for(int z=0; z < coolVec.size(); z++){
		binDensity->FillComponent(0, 0.0);

		diffusionCoeff = coolVec[z][0];
		Cd = 1.0 / diffusionCoeff;
		viterMax = coolVec[z][1];
		currAm = coolVec[z][2];
		percentStrain = coolVec[z][3];
		printStep = (int)coolVec[z][4];
		viscosity = Cd / dt;

		s = (100 / (EdgeLength*percentStrain))*log(2.0);
		bd->updateProperty(OPSBody::aM, currAm);
		bd->updateProperty( OPSBody::sv, s );

		bool checkConsistency = false;
		if (checkConsistency) {
			std::cout << "Checking consistency......" << std::endl;
			bk->updateParallelKick();
			bd->checkConsistency(true);
		}

		std::cout << "Viscosity Input Parameters:" << std::endl
				<< " Cd = " << Cd << std::endl
				<< "  D = " << diffusionCoeff << std::endl
				<< " dt = " << dt << std::endl;

		//***************************  INNER SOLUTION LOOP ***************************//

		for (int viter = 0; viter < viterMax; viter++) {

			cout << endl
					<< "VISCOUS ITERATION: " << viter + stepCount
					<< "\t viscosity = " << vr->viscosity()
					<< endl
					<< endl;

			bk->updateParallelKick();
			//std::cout<<"Average kick norm: "<< bk.getKickStats()
			//<<std::endl;

			bool printBeforeSolve = false;
			if (printBeforeSolve) {
				//We will print only after every currPrintStep iterations
				if (viter % printStep == 0) {
					sstm << fname << "-initial-" << nameSuffix;
					rName = sstm.str();
					bd->printParaview(rName.c_str());
					sstm.str("");
					sstm.clear();
				}
			}

			solver.solve(&model);

			std::cout << "VISCOSITY: " << std::endl
					<< "          velocity = " << vr->velocity() << std::endl
					<< " updated viscosity = " << vr->viscosity() << std::endl
					<< std::endl;

			//Fill-in Matrix for new state for Kabsch algorithm
			for (int col = 0; col < nodes.size(); col++) {
				Vector3D coord = nodes[col]->deformedPosition();
				for (int row = 0; row < 6; row++) {
					current(row, col) = coord(row);
				}
			}

			//********** Print relaxed configuration ************//
			Affine6d A;
			Matrix6Xd newCurr(6, nodes.size());
			A = Find6DAffineTransform(current, initial);
			for (int col = 0; col < current.cols(); col++) {
				newCurr.col(col) = A.linear()*current.col(col)
														+ A.translation();
			}

			//Update the Kabsch transformed positions as new current
			//configuration in the node container
			for (int i = 0; i < nodes.size(); i++) {
				Vector6D kabschPoint;
				for (int j = 0; j < 6; j++) {
					kabschPoint(j) = newCurr(j, i);
				}
				nodes[i]->setDeformedDOFs(kabschPoint);
			}
			bd->updatePolyDataAndKdTree();
			bd->updateNeighbors();

			//We will print only after every currPrintStep iterations
			if (viter % printStep == 0) {
				paraviewStep++;
				sstm << fname << "-relaxed-" << nameSuffix++;
				rName = sstm.str();
				bd->printParaview(rName.c_str());
				sstm.str("");
				sstm.clear();
			}

			int paraviewStepPrint;
			paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;

			myfile << z << "\t"
					<< paraviewStepPrint <<"\t"
					<< diffusionCoeff <<"\t"
					<< currAm << "\t"
					<< bd->getAsphericity() << "\t"
					//<< bd->getLoopAsphericity() << "\t"
					<< bd->getAverageRadius() << "\t"
					<< bd->getMorseEnergy() << "\t"
					<< bd->getPlanarityEnergy() << "\t"
					<< bd->getNormalityEnergy() << "\t"
					<< bd->getCircularityEnergy() << "\t"
					<< bk->energy() << "\t"
					<< vr->energy() << "\t"
					<< solver.function()
					<< endl;
			//********** Find bins for each particle ************//
			putParticlesInBins(cellLimits, newCurr, nodes.size(), binDensity, viterMax);

			// step forward in "time", relaxing viscous energy & forces
			vr->step();
		}
	}
	myfile.close();
	t2 = clock();
	float diff((float)t2 - (float)t1);
	std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
			<< " seconds" << std::endl;

	//Release the dynamically allocated memory
	delete bd;
	delete bk;
	delete vr;
	for(vector<OPSNode*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i){
		delete *i;
	}
	nodes.clear();
}
