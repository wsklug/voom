#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Contact.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include "HelperFunctions.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
	clock_t t1, t2;
	t1 = clock();
	if (argc < 2) {
		cout << "Usage: indent modelName"
			<< endl;
		return(0);
	}
	bool verbose = true;
	string modelName = argv[1];

	double gamma_inp = 0.0;
	double indent_inp = 0.0;
	double viscosity_inp = 0.0;
	double friction_inp = 0.0;
	double step_inp = 0.0;
	bool remesh = true;
	bool unload = false;
	double Ravg = 0;
	double Zmin;
	double Zmax;
	double afmR;
	tvmet::Vector<double, 3> xc(0.0);
	double Z_glass;
	double dZ;
	double ARtol = 1.05;
	double cleanTol = 0.01;
	double capsoSearchRadFactor = 1.5;
    int loopSurfSubDiv = 3;

	//Read epsilon and percentStrain from input file.
	string temp;
	std::ifstream miscInpFile("miscInp.dat");
	assert(miscInpFile);
	miscInpFile >> temp >> gamma_inp
		>> temp >> indent_inp
		>> temp >> step_inp
		>> temp >> friction_inp
		>> temp >> cleanTol
		>> temp >> capsoSearchRadFactor
		>> temp >> loopSurfSubDiv
		>> temp >> unload;

	miscInpFile.close();

	if (gamma_inp <= 0.0) {
		std::cout << "gamma = " << gamma_inp << " but should be positive."
			<< std::endl;
		return 0;
	}

	vtkPolyDataReader * reader = vtkPolyDataReader::New();
	vtkSmartPointer<vtkPolyData> mesh;
	string inputFileName;
	vtkSmartPointer<vtkDataArray> displacements;

	inputFileName = modelName + ".vtk";
	reader->SetFileName(inputFileName.c_str());
	//We will use this object, shortly, to ensure consistent triangle orientations
	vtkPolyDataNormals * normals = vtkPolyDataNormals::New();

	//We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput()
	vtkSmartPointer<vtkPolyData> ds = reader->GetOutput();
	reader->Update();
	normals->SetInputConnection(reader->GetOutputPort());
	// send through normals filter to ensure that triangle orientations
	// are consistent
    normals->ComputeCellNormalsOn();
	normals->ConsistencyOn();
	normals->SplittingOff();
	normals->AutoOrientNormalsOn();

	mesh = normals->GetOutput();
	normals->Update();
	std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
		<< std::endl;

	// create vector of nodes
	int dof = 0;
	std::vector< NodeBase* > nodes;
	std::vector< DeformationNode<3>* > defNodes;

	// read in points
	for (int a = 0; a < mesh->GetNumberOfPoints(); a++) {
		int id = a;
		DeformationNode<3>::Point x;
		mesh->GetPoint(a, &(x[0]));
		Ravg += tvmet::norm2(x);
		NodeBase::DofIndexMap idx(3);
		for (int j = 0; j < 3; j++) idx[j] = dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
		nodes.push_back(n);
		defNodes.push_back(n);
	}

	assert(nodes.size() != 0);
	Ravg /= nodes.size();
	cout << "Number of nodes: " << nodes.size() << endl
		<< "Ravg = " << Ravg << endl;

	// read in triangle connectivities
	vector< tvmet::Vector<int, 3> > connectivities;
	tvmet::Vector<int, 3> c;
	int ntri = mesh->GetNumberOfCells();
	connectivities.reserve(ntri);
	if (verbose) cout << "Number of triangles: " << ntri << endl;
	for (int i = 0; i < ntri; i++) {
		assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
		for (int a = 0; a < 3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
		connectivities.push_back(c);
	}

	std::vector<double> lengthStat;
	double EdgeLength;
	double stdDevEdgeLen;
	double C0 = 0.0;

	// Calculate side lengths average and std dev of the 
	//equilateral triangles
	lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
	EdgeLength = lengthStat[0];
	stdDevEdgeLen = lengthStat[1];
	std::cout << "Before any relaxation :" << endl
		<< "   Average triangle edge length = " << std::setprecision(10)
		<< EdgeLength << endl
		<< "   Standard deviation = " << std::setprecision(10)
		<< stdDevEdgeLen << endl;
	std::cout.precision(6);

	//Rescale the capsid such that triangle edge-lengths are unity
	for (int i = 0; i < defNodes.size(); i++) {
		DeformationNode<3>::Point x;
		x = defNodes[i]->point();
		x *= 1.0 / EdgeLength;
		defNodes[i]->setPoint(x);
		defNodes[i]->setPosition(x);
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

	//Set the AFM indenter radius
	afmR = Ravg;

	//Material properties
	double gamma = gamma_inp;
	double Y = std::sqrt(gamma) / Ravg; //2D Young's modulus
	double nu = 1.0 / 3.0;
	double KC = Ravg / std::sqrt(gamma); //Bending modulus
	double KG = -2 * (1 - nu)*KC; // Gaussian modulus
	C0 = 0.0;
	int quadOrder = 2;

	int m = 5;
	int maxIter = 1e5;
	double factr = 1.0e+1;
	double pgtol = 1.0e-7;
	int iprint = 1000;

	std::stringstream sstm;
	string fname = modelName;
	string rName;

	typedef FVK MaterialType;
	typedef LoopShellBody<MaterialType> LSB;
	typedef LoopShell<MaterialType> LS;
	MaterialType bending(KC, KG, C0, Y, 0.3);

	//********************************** Actual Indentation **********************//

	LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
	bd->setOutput(paraview);

	// create Model
	Model::BodyContainer bdc;
	bdc.push_back(bd);

	Model model(bdc, nodes);

	Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, 1e5);
	std::cout << "Relaxing shape for gamma = " << gamma << std::endl;

	// relax initial shape;
	solver.solve(&model);

	std::cout << "Shape relaxed." << std::endl
		<< "Energy = " << solver.function() << std::endl;

	//Calculate centre of sphere as average of position vectors of all nodes.
	tvmet::Vector<double, 3> Xavg(0.0);
	for (int i = 0; i < defNodes.size(); i++) {
		Xavg += defNodes[i]->point();
	}
	Xavg /= defNodes.size();

	//We will calculate radius using the quadrature points
	vtkSmartPointer<vtkPolyData> lssPd = 
        bd->getLoopShellSurfPoints(cleanTol, loopSurfSubDiv);
	std::vector<double> radialStats = getRadialStats(lssPd, Xavg);
	Ravg = radialStats[0];
	std::cout << "Radius of capsid after relaxation = " << Ravg << endl;
	double asphericity = radialStats[1];
	double capsomerSearchRad = radialStats[2];

	std::cout << "Effective 2D Young's modulus = " << Y << endl
		<< "FVK number = " << gamma_inp << endl
		<< "Asphericity = " << asphericity << endl;

	// find top and bottom of capsid
	Zmin = std::numeric_limits<double>::max();
	Zmax = -std::numeric_limits<double>::max();

	for (int a = 0; a < defNodes.size(); a++) {
		double Z = defNodes[a]->getPoint(2);
		Zmin = std::min(Zmin, Z);
		Zmax = std::max(Zmax, Z);
	}

	std::cout << "Zmax = " << Zmax << std::endl
		<< "Zmin = " << Zmin << std::endl;

	xc = 0.0, 0.0, Zmax + afmR;

	std::cout << "AFM Indenter radius =" << afmR << std::endl
		<< "AFM Indenter center = (" << xc[0] << ","
		<< xc[1] << "," << xc[2] << ")" << std::endl;

	std::cout << "Compressing capsid." << std::endl;

	double friction = friction_inp;
	double k_AL = 1.0e2;
	RigidHemisphereAL * afm
		= new RigidHemisphereAL(defNodes, k_AL, afmR, xc, friction);
	afm->updateContact();

	bd->pushBack(afm);
	std::cout << "Added afm to body." << std::endl;

	bool up = true;
	Z_glass = Zmin;

	RigidPlateAL* glass = new RigidPlateAL(defNodes, k_AL, Z_glass, up, friction);
	bd->pushBack(glass);

	const double originalHeight = Zmax - Zmin;
	double Zbegin = Z_glass;
	double Zend = Zmin + 1.0*Ravg;
	if (indent_inp > 0.0) {
		Zend = Zbegin + indent_inp*Ravg;
	}

	dZ = (Zend - Zbegin) / 100;
	if (step_inp > 0.0) {
		dZ = step_inp*Ravg;
	}

	double bdEnergy = bd->energy();

	string fzName = modelName + ".fz";
	ofstream FvsZ(fzName.c_str());
	FvsZ << "#Step\tIndentation\tGlass_Fz\tAFM_Fz"
		<< "\tLSBStrainEnergy\tPlateEnergy\tAFMEnergy\tSolverFunction"
		<< std::endl;

	blitz::Array<double, 1> x_prev(model.dof());
	blitz::Array<double, 1> u_prev(model.dof());

	model.getField(solver);
	for (int i = 0; i < model.dof(); i++) x_prev(i) = solver.field(i);
	u_prev = 0.0;

	//Store all the file names in this vector
	std::vector<string> allStepFiles;

	// %%%%%%%%%%%%%%%%%%%%%%
	// Begin indentation loop
	// %%%%%%%%%%%%%%%%%%%%%%

	int step = 0;
	for (double Z = Zbegin; ; Z += dZ, step++) {

		// initial guess
		if (step == 0) {
			// shift capsid up by dZ/2 as an initial guess
			for (int a = 0; a < defNodes.size(); a++) {
				defNodes[a]->addPoint(2, 0.5*dZ);
			}
			model.getField(solver);
		}
		else if (std::abs(dZ) > 0) { // 
									 // add scaled version of previous displacement as an initial
									 // guess, but only in loading direction
			for (int i = 0; i < model.dof(); i++) {
				solver.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
			}
			model.putField(solver);
		}

		// move glass up by dZ
		glass->setZ(Z);

		std::cout << std::endl
			<< "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
			<< std::endl
			<< " Z = " << Z
			<< " zeta = " << Z - Zbegin
			<< std::endl;

		// update contact
		glass->updateContact();
		afm->updateContact();
		/*
		bool checkConsistency = true;
		if (checkConsistency) {
			std::cout << "Checking consistency......" << std::endl;
			model.checkConsistency(true, false);
		}
		*/
		model.computeAndAssemble(solver, false, true, false);
		bd->printParaview("contact");

		solver.solve(&model);
		bdEnergy = bd->energy();

		if (verbose) {
			std::cout << "ENERGY:" << std::endl
				<< "   body energy = " << bdEnergy << std::endl
				<< "  total energy = " << solver.function() << std::endl
				<< std::endl;
		}

		// update contact
		if (verbose) {
			std::cout << "CONTACT:" << std::endl;
			std::cout << "       top active  = " << afm->active() << std::endl
				<< "    bottom active  = " << glass->active() << std::endl
				<< "        top force  = " << afm->FZ() << std::endl
				<< "     bottom force  = " << glass->FZ() << std::endl
				<< "  top penetration  = " << afm->penetration() << std::endl
				<< "bottom penetration = " << glass->penetration() << std::endl;
		}


		double height = Zmax - Z;
		// add up forces on top and bottom
		double Ztop = Zmax;
		double Zbot = Zmax - height;

		if (verbose) {
			std::cout << "Contact energy converged." << std::endl
				<< std::endl
				<< "        height = " << height << std::endl
				<< "   indentation = " << originalHeight - height << std::endl
				<< std::endl;
		}

		// Compute displacement
		for (int i = 0; i < model.dof(); i++) {
			u_prev(i) = (solver.field(i) - x_prev(i)) / std::abs(dZ);
		}

		// If displacement was bigger than dZ, some buckling event must
		// have occurred and it's better not to make a continuation
		// attempt.
		if (max(abs(u_prev)) > 1.0) u_prev = 0.0;//????

												 // Save current (successful) state as previous
		for (int i = 0; i < model.dof(); i++) x_prev(i) = solver.field(i);

		if (unload && Z >= Zend) { // reached max indentation, now unload by reversing dZ
			dZ = -dZ;
		}

		//*********** BEGIN PRINTING OUTPUT (and log) FILES ***********//

		FvsZ << step
			<< "\t" << (originalHeight - height) / Ravg
			<< "\t" << glass->FZ() / std::sqrt(Y*KC)
			<< "\t" << afm->FZ() / std::sqrt(Y*KC)
			<< "\t" << bd->totalStrainEnergy() / (Ravg*std::sqrt(Y*KC))
			<< "\t" << glass->energy() / (Ravg*std::sqrt(Y*KC))
			<< "\t" << afm->energy() / (Ravg*std::sqrt(Y*KC))
			<< "\t" << solver.function() / (Ravg*std::sqrt(Y*KC))
			<< std::endl;

		sstm << modelName << "-step-" << step;
		rName = sstm.str();
		bd->printParaview(rName.c_str());
		//We will append Caspsomer POINT_DATA to the vtk output file
		//printed by printParaview(), if such a file exists
		sstm << ".vtk";
		rName = sstm.str();
		sstm.str("");
		sstm.clear();
		allStepFiles.push_back(rName);

		//Now we will print the LoopShellSurface
		//We will calculate radius using the quadrature points
		lssPd = bd->getLoopShellSurfPoints(cleanTol, loopSurfSubDiv);
		radialStats = getRadialStats(lssPd, Xavg);
		capsomerSearchRad = radialStats[2];
		sstm << modelName << "-LoopShellSurf-"
			<< step << ".vtk";
		rName = sstm.str();
		meshSphericalPointCloud(lssPd, capsoSearchRadFactor*capsomerSearchRad,
			rName);
		std::cout << "\tCapsomer Search Radius = " << capsomerSearchRad
			<< std::endl;
		sstm.str("");
		sstm.clear();
		//************* END PRINTING OUTPUT FILES **************//

		// check if we are done
		if (unload && Z + dZ < Zbegin)
			break;
		else if (!unload && Z + dZ > Zend)
			break;

	}// Indentation Loop Ends
	insertValenceInVtk(allStepFiles);
	FvsZ.close();

	std::cout << "Indentation complete." << std::endl;

	t2 = clock();
	float diff((float)t2 - (float)t1);
	std::cout << "Total execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
	return 0;

}

///////////////////////// END OF MAIN FUNCTION //////////////////////////////////
/////////////////////////                     //////////////////////////////////
