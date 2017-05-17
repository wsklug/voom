#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Quadrature.h"
#include "TriangleQuadrature.h"

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>

#include "Morse.h"
#include "SpringPotential.h"
#include "PotentialBody.h"
#include "HelperFunctions.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
	clock_t t1, t2, t3;
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

	//For Morse material 
	double epsilon;
	double Rshift;
	double percentStrain = 10;
	bool harmonicRelaxNeeded;
	double searchRadFactor = 1.2;
	double capsoSearchRadFactor = 1.2;
	bool projectOnSphere = true;
	double ARtol = 1.5;
	bool remesh = true;
	bool useOldRemeshTechnicque = false;
	double cleanTol = 0.0;


	//Read epsilon and percentStrain from input file. percentStrain is
	//calculated so as to set the inflection point of Morse potential
	//at a fixed distance relative to the equilibrium separation
	//e.g. 1.1*R_eq, 1.5*R_eq etc.
	std::ifstream miscInpFile("miscInp.dat");
	assert(miscInpFile);
	string temp;
	miscInpFile >> temp >> epsilon
		>> temp >> percentStrain
		>> temp >> ARtol
		>> temp >> searchRadFactor
		>> temp >> cleanTol
		>> temp >> capsoSearchRadFactor
		>> temp >> harmonicRelaxNeeded
		>> temp >> projectOnSphere
		>> temp >> remesh
		>> temp >> useOldRemeshTechnicque;

	miscInpFile.close();

	vtkSmartPointer<vtkDataSetReader> reader =
		vtkSmartPointer<vtkDataSetReader>::New();

	vtkSmartPointer<vtkPolyData> mesh;

	std::stringstream sstm;


	reader->SetFileName(inputFileName.c_str());
	//We will use this object, shortly, to ensure consistent triangle
	//orientations
	vtkSmartPointer<vtkPolyDataNormals> normals =
		vtkSmartPointer<vtkPolyDataNormals>::New();

	//We have to pass a vtkPolyData to vtkPolyDataNormals::SetInputData()
	//If our input vtk file has vtkUnstructuredGridData instead of
	//vtkPolyData then we need to convert it using vtkGeometryFilter
	reader->Update();
	vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();
	normals->SetInputConnection(reader->GetOutputPort());
	// send through normals filter to ensure that triangle orientations
	// are consistent 
	normals->ConsistencyOn();
	normals->SplittingOff();
	normals->AutoOrientNormalsOn();
	mesh = normals->GetOutput();
	normals->Update();
	
	std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
		<< std::endl;

	vtkSmartPointer<vtkDataArray> displacements;

	// create vector of nodes

	int dof = 0;
	std::vector< NodeBase* > nodes;
	std::vector< DeformationNode<3>* > defNodes;
	double Ravg = 0.0;

	// read in points
	for (int a = 0; a < mesh->GetNumberOfPoints(); a++) {

		int id = a;
		NodeBase::DofIndexMap idx(3);

		for (int j = 0; j < 3; j++) idx[j] = dof++;

		DeformationNode<3>::Point X;
		DeformationNode<3>* n;

		mesh->GetPoint(a, &(X[0]));
		n = new DeformationNode<3>(id, idx, X);
		Ravg += tvmet::norm2(X);

		nodes.push_back(n);
		defNodes.push_back(n);
	}
	assert(nodes.size() != 0);
	Ravg /= nodes.size();
	cout << "Number of nodes: " << nodes.size() << endl
		<< "Initial radius: " << Ravg << endl;

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
	std::cout << "Before any relaxation :" << endl
		<< "   Average triangle edge length = "	<< EdgeLength << endl
		<< "   Standard deviation = " << stdDevEdgeLen << endl;
	// Rescale size of the capsid by the average equilateral edge length
	for (int i = 0; i < defNodes.size(); i++) {
		DeformationNode<3>::Point X;
		X = defNodes[i]->position();
		X *= 1.0 / EdgeLength;
		defNodes[i]->setPoint(X);
		defNodes[i]->setPosition(X);
	}

	//Recalculate edge lengths and capsid radius
	lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
	EdgeLength = lengthStat[0];

	Ravg = 0.0;
	for (int i = 0; i < defNodes.size(); i++) {
		DeformationNode<3>::Point x;
		x = defNodes[i]->point();
		double tempRadius = tvmet::norm2(x);
		Ravg += tempRadius;
	}
	Ravg /= defNodes.size();

	std::cout << "Radius of capsid after rescaling = " << Ravg << endl;

	if (projectOnSphere) {

		//Project points to a sphere of radius Ravg
		for (int i = 0; i < defNodes.size(); i++) {
			DeformationNode<3>::Point X;
			X = defNodes[i]->position();
			X *= Ravg / (tvmet::norm2(X));
			defNodes[i]->setPoint(X);
			defNodes[i]->setPosition(X);
		}

		//Recalculate edge lengths and capsid radius
		lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
		EdgeLength = lengthStat[0];

		Ravg = 0.0;
		for (int i = 0; i < defNodes.size(); i++) {
			DeformationNode<3>::Point x;
			x = defNodes[i]->point();
			double tempRadius = tvmet::norm2(x);
			Ravg += tempRadius;
		}
		Ravg /= defNodes.size();

		std::cout << "Radius of capsid after projecting on a sphere = " << Ravg << endl;
	}
	//******************* READ FVK DATA FROM FILE **********************

	// We want variable number of FVK increments in different ranges. So
	// we will read FVK values from a file instead of generating it by
	// code  
	std::ifstream fvkFile("fvkSteps.dat");
	assert(fvkFile);
	std::vector<vector<double> > gammaVec;
	double currFVK, currPrintFlag;
	while (fvkFile >> currFVK >> currPrintFlag) {
		std::vector<double> currLine;
		currLine.push_back(currFVK);
		currLine.push_back(currPrintFlag);
		gammaVec.push_back(currLine);
	}
	fvkFile.close();

	string fname = inputFileName.substr(0, inputFileName.find("."));
	string iName;
	string rName;
	string actualFile;

	ofstream myfile;
	std::string dataOutputFile = "asphVsFVKMorse.dat";

	myfile.open(dataOutputFile.c_str());
	myfile << "#Step" << "\t" << "Ravg" << "\t"
		<< "Y" << "\t" << "FVKin" << "\t"
		<< "asphericity" << "\t" << setw(8) << "LSBStrainEnergy" << "\t"
		<< "MorseEnergy" << "\t"
		<< "TotalFunctional" << "\t"
		<< "RemeshOccured" << "\t"
		<< "LargestARtol" << "\t"
		<< "EnergyBeforeRemesh" << "\t"
		<< "EnergyAfterRemesh" << std::endl;
	myfile << showpoint;


	//Parameters for the l-BFGS solver
	int m = 7;

	//int maxIter = 1e6;
	int maxIter1 = 1e5;
	int maxIter2 = 1e5;

	double factr = 10.0;
	double pgtol = 1e-7;
	int iprint = 1000;
	Lbfgsb solver1(3 * nodes.size(), m, factr, pgtol, iprint, maxIter1, true);

	typedef FVK MaterialType;
	typedef LoopShellBody<MaterialType> LSB;
	typedef LoopShell<MaterialType> LS;

	//****************  Protein body parameters *******************

	Rshift = EdgeLength;
	double sigma = (100 / (Rshift*percentStrain))*log(2.0);

	double PotentialSearchRF;
	double springConstant = 2 * sigma*sigma*epsilon;
	double gamma = gammaVec[0][0];
	double Y = 2.0 / sqrt(3)*springConstant; // Young's modulus
	double nu = 1.0 / 3.0;
	double KC = Y*Ravg*Ravg / gamma;
	double KG = -2 * (1 - nu)*KC; // Gaussian modulus
	double C0 = 0.0;
	int quadOrder = 2;

	if (harmonicRelaxNeeded) {

		//****** Relax the initial mesh using harmonic potential ****** //

		MaterialType bending(KC, KG, C0, 0.0, 0.0);
		LSB * bd1 = new LSB(bending, connectivities, nodes, quadOrder);
		bd1->setOutput(paraview);
		SpringPotential SpringMat(springConstant, Rshift);
		PotentialSearchRF = 1.2*Rshift;
		PotentialBody * SpringBody = new
			PotentialBody(&SpringMat, defNodes, PotentialSearchRF);

		//Create Model
		Model::BodyContainer bdc1;
		bdc1.push_back(SpringBody);
		bdc1.push_back(bd1);
		Model model1(bdc1, nodes);

		std::cout << "Spring constant: " << springConstant << endl;
		std::cout << "Relaxing the mesh using harmonic potential..." << endl;
		solver1.solve(&model1);

		std::cout << "Harmonic potential relaxation complete." << endl;

		//Print to VTK file
		/*
		sstm << fname << "-relaxed-harmonic";
		rName = sstm.str();
		bd1->printParaview(rName);
		sstm.str("");
		sstm.clear();
		*/
		//Release allocated memory
		delete bd1;
		delete SpringBody;

		//Recalculate edge lengths and dependent quantities
		lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
		EdgeLength = lengthStat[0];
		Rshift = EdgeLength;
		sigma = (100 / (Rshift*percentStrain))*log(2.0);
		springConstant = 2 * sigma*sigma*epsilon;
		std::cout << "After relaxing with harmonic potential: " << endl
			<< "   Average triangle edge length = "	<< EdgeLength << endl
			<< "   Standard deviation = " << stdDevEdgeLen << endl;		
	}

	//***************************  SOLUTION LOOP ***************************

	//Loop over all values of gamma and relax the shapes to get
	//asphericity
	std::cout << "Morse potential parameters:" << endl
		<< "sigma = " << sigma << " epsilon = " << epsilon
		<< " Rshift = " << Rshift << endl;

	PotentialSearchRF = searchRadFactor*Rshift;
	std::cout << "PotentialSearchRF = " << PotentialSearchRF
		<< std::endl;
	//Morse parameters do not change with FvK
	Morse Mat(epsilon, sigma, Rshift);
	PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);
	//Assign the MorseBond structure from the initial neighbor information
	vtkSmartPointer<vtkCellArray> bonds = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType bond[2] = { 0, 0 };
	for (int i = 0; i < defNodes.size(); i++) {
		tvmet::Vector<double, 3> centerNode = defNodes[i]->point();
		for (int a = i + 1; a < defNodes.size(); a++) {
			tvmet::Vector<double, 3> currNode = defNodes[a]->point();
			double r = tvmet::norm2(centerNode - currNode);
			if (r <= PotentialSearchRF) {
				bond[0] = i;
				bond[1] = a;
				bonds->InsertNextCell(2, bond);
			}
		}
	}
	for (int q = 0; q < gammaVec.size(); q++) {

		gamma = gammaVec[q][0];
		std::cout << "########### Gamma = " << gamma << "###########"
			<< std::endl;
		double currPrintFlag = gammaVec[q][1];

		// Bending body parameters
		Y = 2.0 / sqrt(3)*springConstant; // Young's modulus
		KC = Y*Ravg*Ravg / gamma; // Bending modulus
		KG = -2 * (1 - nu)*KC; // Gaussian modulus

		//The Bodies
		MaterialType bending(KC, KG, C0, 0.0, 0.0);
		LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
		bd->setOutput(paraview);

		//Create Model
		Model::BodyContainer bdc;
		bdc.push_back(PrBody);
		bdc.push_back(bd);
		Model model(bdc, nodes);

		double energy = 0;

		solver1.solve(&model);
		energy = solver1.function();

		// REMESHING
		uint elementsChanged = 0;
		bool remeshEventFlag = false;
		int numAttempts = 0;
		double energy_prev = energy;
		std::vector<double> meshQuality = { 0.0,0.0,0.0 };
		if (remesh) {
			meshQuality = getMeshQualityInfo(defNodes, connectivities);
			if (useOldRemeshTechnicque) {
				elementsChanged = bd->Remesh(ARtol, bending, quadOrder);
				//Print out the number of elements that changed due to remeshing
				if (elementsChanged > 0) {
					remeshEventFlag = true;
					std::cout << "Number of elements that changed after remeshing = "
						<< elementsChanged << "." << std::endl;
					//Relax again after remeshing
					solver1.solve(&model);
					meshQuality = getMeshQualityInfo(defNodes, connectivities);
					std::cout << "\tLargest aspect ratio after re-solving = "
						<< meshQuality[2] << std::endl;
					std::cout << "Energy before remeshing = " << energy << std::endl;
					energy = solver1.function();
					std::cout << "Energy after remeshing  = " << energy << std::endl;
				}
			}
			else {				
				while (meshQuality[2] > ARtol && numAttempts < 1) {
					remeshEventFlag = true;
					std::cout << "***** Skewed triangles detected in the mesh. *****" << std::endl
						<< "\tLargest aspect ratio before remeshing = "
						<< meshQuality[2] << std::endl;
					std::cout << "Remeshing Attempt # " << numAttempts << std::endl;
					//Try remeshing
					connectivities = delaunay3DSurf(defNodes);
					meshQuality = getMeshQualityInfo(defNodes, connectivities);
					std::cout << "\tLargest aspect ratio after remeshing = "
						<< meshQuality[2] << std::endl;
					bd->CreateLoopFE(connectivities, bending, quadOrder, true);
					solver1.solve(&model);
					std::cout << "Energy before remeshing = " << energy << std::endl;
					energy = solver1.function();
					std::cout << "Energy after remeshing  = " << energy << std::endl;
					meshQuality = getMeshQualityInfo(defNodes, connectivities);
					std::cout << "\tLargest aspect ratio after re-solving = "
						<< meshQuality[2] << std::endl;					
					numAttempts++;
				}
			}
		}// Remeshing ends here

		//Calculate maximum principal strains in all elements
		//bd->calcMaxPrincipalStrains();

		std::cout << "Shape relaxed." << std::endl
			<< "Energy = " << energy << std::endl;		

		//Re-calculate triangle edge length mean and deviation
		lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
		EdgeLength = lengthStat[0];
		stdDevEdgeLen = lengthStat[1];
		std::cout << "After relaxing with Morse potential: " << endl
			<< "   Average triangle edge length = " << EdgeLength << endl
			<< "   Standard deviation = " << stdDevEdgeLen << endl;

		//Calculate centre of sphere as average of position vectors of all nodes.
		tvmet::Vector<double, 3> Xavg(0.0);
		for (int i = 0; i < defNodes.size(); i++) {
			Xavg += defNodes[i]->point();
		}
		Xavg /= defNodes.size();

		//We will calculate radius using the quadrature points
		vtkSmartPointer<vtkPolyData> lssPd = bd->getLoopShellSurfPoints(cleanTol);
		std::vector<double> radialStats = getRadialStats(lssPd, Xavg);
		Ravg = radialStats[0];
		std::cout << "Radius of capsid after relaxation = " << Ravg << endl;
		double asphericity = radialStats[1];
		double capsomerSearchRad = radialStats[2];		
		//Selectively print the relaxed shapes
		if (currPrintFlag) {
			sstm << fname << "-relaxed-" << q;
			rName = sstm.str();
			bd->printParaview(rName);
			sstm.str("");
			sstm.clear();
			//Now we will print the LoopShellSurface
			sstm << fname << "-LoopShellSurf-"
				<< q << ".vtk";
			rName = sstm.str();			
			meshSphericalPointCloud(lssPd, capsoSearchRadFactor*capsomerSearchRad,
				rName);
			std::cout << "\tCapsomer Search Radius = " << capsomerSearchRad
				<< std::endl;
			sstm.str("");
			sstm.clear();
		}

		myfile << setprecision(16) << q << "\t" << Ravg << "\t" << Y
			<< "\t" << gamma << "\t" << asphericity << "\t" 
			<< bd->totalStrainEnergy() << "\t" << PrBody->energy() << "\t" 
			<< energy << "\t" << remeshEventFlag << "\t" << meshQuality[2] 
			<< "\t" << energy_prev << "\t" << "\t" << energy << std::endl;		

		//Release the dynamically allocated memory
		delete bd;
	}
	delete PrBody;
	std::vector< NodeBase* >().swap(nodes);
	std::vector< DeformationNode<3>* >().swap(defNodes);
	myfile.close();
	t2 = clock();
	float diff((float)t2 - (float)t1);
	std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

	// Post-processing: Manipulating VTK files
	std::vector<std::string> allVTKFiles;
	int numFiles = gammaVec.size();
	allVTKFiles.reserve(numFiles);
	for (int fileNum = 0; fileNum < numFiles; fileNum++) {
		sstm << fname << "-relaxed-" << fileNum << ".vtk";
		std::string tempString = sstm.str();
		allVTKFiles.push_back(tempString);
		sstm.str("");
		sstm.clear();
	}

	//insertValenceInVtk(allVTKFiles);
	//writeEdgeStrainVtk(allVTKFiles, Rshift, percentStrain);
	plotMorseBonds(allVTKFiles, fname, epsilon, Rshift, sigma, bonds);
	t3 = clock();
	diff = ((float)t3 - (float)t2);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

}
