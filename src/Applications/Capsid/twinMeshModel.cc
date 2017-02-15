#include <tvmet/Vector.h>
#include <iomanip>
#include <limits>
#include "Node.h"
#include "FVK.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Quadrature.h"
#include "LoopShellBody.h"
#include "TriangleQuadrature.h"

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>
#include <vtkSphereSource.h>
#include <vtkCellArray.h>
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>

#include "Morse.h"
#include "SpringPotential.h"
#include "PotentialBody.h"
#include "BrownianKick.h"
#include "ViscousRegularizer.h"
#include "RadialSpring.h"
#include "ViscosityBody.h"

#include "HelperFunctions.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;


typedef DeformationNode<3> Node; // nickname for mechanics nodes
typedef std::vector< Node* > NodeContainer;

int main(int argc, char* argv[]) {

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
	double percentStrain;
	double pressureFactor;

	//Numerical viscosity input parameter
	int viterMax;

	double Rshift;
	bool rescale;
	bool areaConstraintOn = false;
	bool pressureConstraintOn = false;
	double dt = 9.76e-4;
	int nameSuffix = 0;
	int lat_res = 100, long_res = 101;
	double scaleC0 = 1.0;
	int numSubDivide = 0.0;

	//Auxiliary input file
	std::ifstream miscInpFile("miscInp.dat");
	assert(miscInpFile);
	string temp;
	miscInpFile >> temp >> epsilon
		>> temp >> pressureFactor
		>> temp >> dt
		>> temp >> pressureConstraintOn
		>> temp >> areaConstraintOn
		>> temp >> rescale
		>> temp >> lat_res
		>> temp >> long_res
		>> temp >> scaleC0
		>> temp >> numSubDivide;

	miscInpFile.close();

	std::stringstream sstm;

	vtkSmartPointer<vtkDataSetReader> reader =
		vtkSmartPointer<vtkDataSetReader>::New();
	reader->SetFileName(inputFileName.c_str());

	//Re-orient normals if needed
	vtkSmartPointer<vtkPolyDataNormals> normals =
		vtkSmartPointer<vtkPolyDataNormals>::New();
	normals->SetInputConnection(reader->GetOutputPort());
	normals->ConsistencyOn();
	normals->SplittingOff();
	normals->AutoOrientNormalsOn();
	normals->Update();

	vtkSmartPointer<vtkPolyData> mesh = normals->GetOutput();
	std::cout << "Mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
		<< std::endl;

	//Create a finer mesh by sub-dividing triangles from the original mesh
	vtkSmartPointer<vtkTriangleFilter> triangles
		= vtkSmartPointer<vtkTriangleFilter>::New();
	triangles->SetInputData(mesh);

	vtkSmartPointer<vtkLinearSubdivisionFilter> linSub
		= vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
	//vtkSmartPointer<vtkLoopSubdivisionFilter> linSub
		//= vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
	linSub->SetInputConnection(triangles->GetOutputPort());
	linSub->SetNumberOfSubdivisions(numSubDivide);
	linSub->Update();

	vtkSmartPointer<vtkPolyData> finerMesh = linSub->GetOutput();

	//Just for testing.. print the finer mesh
	vtkSmartPointer<vtkPolyDataWriter> writer
		= vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(finerMesh);
	writer->SetFileName("FinerMesh.vtk");
	writer->Write();

	//Create vector of nodes
	int dof = 0;
	std::vector< DeformationNode<3>* > defNodes;
	std::vector< NodeBase* > allNodes;
	std::vector< DeformationNode<3>* > allDefNodes;
	double Ravg = 0.0;

	// read in points
	for (int a = 0; a < mesh->GetNumberOfPoints(); a++) {
		int id = a;
		DeformationNode<3>::Point x;
		mesh->GetPoint(a, &(x[0]));
		Ravg += tvmet::norm2(x);
		NodeBase::DofIndexMap idx(3);
		for (int j = 0; j < 3; j++) idx[j] = dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
		defNodes.push_back(n);
		allNodes.push_back(n);
		allDefNodes.push_back(n);
	}
	std::cout << "Number of points in finer mesh:" << finerMesh->GetNumberOfPoints() << std::endl;
	std::cout << "Number of points in coarser mesh:" << mesh->GetNumberOfPoints() << std::endl;

	for (int a = mesh->GetNumberOfPoints(); a < finerMesh->GetNumberOfPoints(); a++) {
		int id = a;
		DeformationNode<3>::Point x;
		finerMesh->GetPoint(a, &(x[0]));
		Ravg += tvmet::norm2(x);
		NodeBase::DofIndexMap idx(3);
		for (int j = 0; j < 3; j++) idx[j] = dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
		allNodes.push_back(n);
		allDefNodes.push_back(n);
	}
	assert(allNodes.size() != 0);
	Ravg /= allNodes.size();
	cout << "Initial radius: " << Ravg << endl;

	// read in finer mesh triangle connectivities
	vector< tvmet::Vector<int, 3> > fineConnectivities;
	int ntri = finerMesh->GetNumberOfCells();
	tvmet::Vector<int, 3> c;
	fineConnectivities.reserve(ntri);
	std::cout << "Number of triangles in finer mesh: " << ntri << endl;

	for (int i = 0; i < ntri; i++) {
		assert(finerMesh->GetCell(i)->GetNumberOfPoints() == 3);
		for (int a = 0; a < 3; a++) c[a] = finerMesh->GetCell(i)->GetPointId(a);
		fineConnectivities.push_back(c);
	}

	// Calculate side lengths average and std dev of the 
	//equilateral triangles
	std::vector<double> lengthStat =
		calcEdgeLenAndStdDev(allDefNodes, fineConnectivities);
	double EdgeLength = lengthStat[0];
	double stdDevEdgeLen = lengthStat[1];
	std::cout << "Before any relaxation :" << endl
		<< "   Average triangle edge length = " << std::setprecision(10)
		<< EdgeLength << endl
		<< "   Standard deviation = " << std::setprecision(10)
		<< stdDevEdgeLen << endl;
	std::cout.precision(6);

	if (rescale) {
		// Rescale size of the capsid by the average equilateral edge length
		for (int i = 0; i < allNodes.size(); i++) {
			DeformationNode<3>::Point x;
			x = allDefNodes[i]->point();
			x *= 1.0 / EdgeLength;
			allDefNodes[i]->setPoint(x);
			allDefNodes[i]->setPosition(x);
		}
		//Recalculate edge lengths and dependent quantities
		lengthStat = calcEdgeLenAndStdDev(allDefNodes, fineConnectivities);
		EdgeLength = lengthStat[0];
		Ravg = 0.0;
		for (int i = 0; i < allDefNodes.size(); i++) {
			DeformationNode<3>::Point x;
			x = allDefNodes[i]->point();
			double tempRadius = tvmet::norm2(x);
			Ravg += tempRadius;
		}
		Ravg /= allDefNodes.size();
		std::cout << "Radius of capsid after rescaling = " << Ravg << endl;
	}

	//Works only for linear subdivision and not for Loop subdivision
	Rshift = EdgeLength*std::pow(2, numSubDivide);

	// Prepare Eigen matrices
	Eigen::Matrix3Xd initial(3, allDefNodes.size()), current(3, allDefNodes.size());

	//Fill in matrix of the initial state for Kabsch algorithm
	for (int col = 0; col < allDefNodes.size(); col++) {
		Vector3D coord = allDefNodes[col]->point();
		for (int row = 0; row < 3; row++) {
			initial(row, col) = coord(row);
		}
	}


	//******************* READ COOLING SCHEDULE from File *******

	std::ifstream coolFile("cooling.dat");
	assert(coolFile);
	std::vector<vector<double> > coolVec;
	double curr_D, currFVK, currPercentStrain, currViterMax, currPrintStep, currRadialSpring;

	std::string headerline;
	std::getline(coolFile, headerline);

	while (coolFile >> curr_D >> currViterMax >> currFVK >>
		currPercentStrain >> currPrintStep >> currRadialSpring) {
		std::vector<double> currLine;
		currLine.push_back(curr_D);
		currLine.push_back(currViterMax);
		currLine.push_back(currFVK);
		currLine.push_back(currPercentStrain);
		currLine.push_back(currPrintStep);
		currLine.push_back(currRadialSpring);
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
		<< "\t" << "Ravg"
		<< "\t" << "asphericity" << "\t" << "FVK" << "\t"
		<< "BendEnergy" << "\t" << "StretchEnergy" << "\t"
		<< "SpringEnergy" << "\t" << "BrownEnergy" << "\t"
		<< "ViscousEnergy" << "\t" << "Total Functional" << "\t"
		<< "MeanSquareDisp"
		<< endl;

	//Parameters for the l-BFGS solver
	int m = 5;
	int maxIter = 1e5;

	double factr = 1.0e+1;
	double pgtol = 1e-7;
	int iprint = 100;
	Lbfgsb solver(3 * allDefNodes.size(), m, factr, pgtol, iprint, maxIter);

	typedef FVK MaterialType;
	typedef LoopShellBody<MaterialType> LSB;
	typedef LoopShell<MaterialType> LS;

	//****************  Protein body parameters ****************//

	double sigma;
	double PotentialSearchRF;
	double springConstant;

	double pressure = 0.0;
	double fracturePressure;

	//The following values have to be choosen as per some
	//guidelines. Read BrownianParameters.pdf to know more

	double gamma;
	double diffusionCoeff = 0;;
	double Cd = 0;
	double viscosity = 0;
	double Y, nu = 1.0 / 3.0;
	double KC = 0, KG = 0, C0 = 0;
	double radialSpringConstant = 0.0;
	int quadOrder = 2;

	double rsEnergy;
	double bkEnergy;
	double vrEnergy;
	double bdEnergy;
	double PrEnergy;
	double energy = 0.0;
	int printStep, stepCount = 0;
	int paraviewStep = -1;

	bool debug = false;

	std::vector< std::string > allVTKFiles;
	std::vector< double > percentStrainData;

	vtkSmartPointer<vtkSphereSource> sp =
		vtkSmartPointer<vtkSphereSource>::New();
	sp->SetRadius(1.0);
	sp->SetThetaResolution(lat_res);
	sp->SetPhiResolution(long_res);
	sp->LatLongTessellationOn();

	vtkSmartPointer<vtkPolyData> pd =
		sp->GetOutput();
	sp->Update();

	std::vector<vector<double> > cellLimits = getSphCellLimits(pd, long_res);

	vtkSmartPointer<vtkDoubleArray> binDensity =
		vtkSmartPointer<vtkDoubleArray>::New();
	binDensity->SetNumberOfComponents(1);
	binDensity->SetNumberOfTuples(pd->GetNumberOfCells());
	binDensity->SetName("Density");

	//***************** MAIN LOOP ****************************//

	for (int q = 0; q < coolVec.size(); q++) {

		binDensity->FillComponent(0, 0.0);

		diffusionCoeff = coolVec[q][0];
		Cd = 1.0 / diffusionCoeff;
		viterMax = coolVec[q][1];
		gamma = coolVec[q][2];
		percentStrain = coolVec[q][3];
		printStep = (int)coolVec[q][4];
		radialSpringConstant = coolVec[q][5];

		viscosity = Cd / dt;

		sigma = (100 / (Rshift*percentStrain))*log(2.0);
		PotentialSearchRF = 1.5*Rshift;
		springConstant = 2 * sigma*sigma*epsilon;

		if (pressureConstraintOn) {
			pressure = 12 * sigma*epsilon
				*(exp(-2 * sigma*Rshift) - exp(-sigma*Rshift)
					+ exp(-1.46410*sigma*Rshift) - exp(-0.7321*sigma*Rshift))
				/ (3 * Ravg*Ravg);

			if (pressure < 0.0) {
				pressure = pressure*(-1);
			}

			fracturePressure = (3.82)*sigma*epsilon / (Rshift*Rshift);
			std::cout << "Fracture Pressure = " << fracturePressure << endl
				<< "Minimum Pressure = " << pressure << endl;
			pressure *= pressureFactor;
			std::cout << "Pressure in use = " << pressure << endl;
		}

		std::cout << "Viscosity Input Parameters:" << std::endl
			<< " Cd = " << Cd << std::endl
			<< "  D = " << diffusionCoeff << std::endl
			<< " dt = " << dt << std::endl;

		Y = 2.0 / sqrt(3)*springConstant; // Young's modulus
		nu = 1.0 / 3.0;
		KC = Y*Ravg*Ravg / gamma;
		KG = -2 * (1 - nu)*KC; // Gaussian modulus
		bool useSpontaneousCurvature = true;

		C0 = useSpontaneousCurvature ? scaleC0*(2 / Ravg) : 0.0;
		std::cout << "Spontaneous Curvature = " << C0 << std::endl;

		//The Bodies
		MaterialType bending(KC, KG, C0, 0.0, 0.0);

		LSB * bd;
		
		std::cout << "Morse potential parameters:" << endl
			<< "sigma = " << sigma << " epsilon = " << epsilon
			<< " Rshift = " << Rshift << endl;

		std::cout << "Pressure = " << pressure << endl
			<< "Capsid radius = " << Ravg << endl;


		//***************************  INNER SOLUTION LOOP ***************************//  

		for (int viter = 0; viter < viterMax; viter++) {
			fineConnectivities = delaunay3DSurf(allDefNodes);
			Morse Mat(epsilon, sigma, Rshift);
			if (areaConstraintOn) {
				std::cout << "********** AREA and PRESSURE CONSTRAINTS ACTIVE  **********"
					<< std::endl;
				bd = new LSB(bending, fineConnectivities, allNodes, quadOrder, pressure,
					0.0, 0.0, 1.0e4, 1.0e6, 1.0e4, multiplier, penalty, noConstraint);
				std::cout << "Prescribed Area = " << bd->prescribedArea() << std::endl;
			}
			else if (pressureConstraintOn && !areaConstraintOn) {
				std::cout << "********** ONLY PRESSURE CONSTRAINT ACTIVE **********" << std::endl;
				bd = new LSB(bending, fineConnectivities, allNodes, quadOrder, pressure,
					0.0, 0.0, 1.0e4, 1.0e6, 1.0e4, multiplier, noConstraint, noConstraint);
			}
			else {
				std::cout << "********** CONSTRAINTS NOT BEING USED **********" << std::endl;
				bd = new LSB(bending, fineConnectivities, allNodes, quadOrder);
			}

			bd->setOutput(paraview);
			PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);
			PrBody->initialNearestNeighbor();
			ViscousRegularizer vr(allNodes, viscosity);
			bd->pushBack(&vr);
			BrownianKick bk(allDefNodes, Cd, diffusionCoeff, dt);
			bd->pushBack(&bk);
			RadialSpring rs(allDefNodes, radialSpringConstant, Ravg);
			bd->pushBack(&rs);

			//Create Model
			Model::BodyContainer bdc;
			bdc.push_back(PrBody);
			bdc.push_back(bd);

			Model model(bdc, allNodes);

			bd->pushBack(&vr);
			bd->pushBack(&bk);
			bd->pushBack(&rs);


			std::cout << std::endl
				<< "VISCOUS ITERATION: " << viter + stepCount
				<< "\t viscosity = " << vr.viscosity()
				<< std::endl
				<< std::endl;

			//bk.updateParallelKick();
			bk.updateProjectedKick();
			std::vector<double> kickStats = bk.getKickStats();
			std::cout << "Average kick norm: " << kickStats[2] << std::endl
				<< "Largest kick norm: " << kickStats[0]
				<< std::endl;

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

			//Fill-in Matrix for new state for Kabsch algorithm
			for (int col = 0; col < allDefNodes.size(); col++) {
				Vector3D coord = allDefNodes[col]->point();
				for (int row = 0; row < 3; row++) {
					current(row, col) = coord(row);
				}
			}

			//We also need to recompute the neighbors for PotentialBody
			PrBody->recomputeNeighbors(PotentialSearchRF);

			rsEnergy = rs.energy();
			bkEnergy = bk.energy();
			vrEnergy = vr.energy();
			bdEnergy = bd->energy() - rsEnergy - vrEnergy - bkEnergy;
			PrEnergy = PrBody->energy();
			energy = solver.function();

			std::cout << "ENERGY:" << std::endl
				<< "viscous energy  = " << vrEnergy << std::endl
				<< "Brownian energy = " << bkEnergy << std::endl
				<< "Spring energy   = " << rsEnergy << std::endl
				<< "protein energy  = " << PrEnergy << std::endl
				<< "bending energy  = " << bdEnergy << std::endl
				<< "  total energy  = " << energy << std::endl
				<< std::endl;
			std::cout << "VISCOSITY: " << std::endl
				<< "          velocity = " << vr.velocity() << std::endl
				<< " updated viscosity = " << vr.viscosity() << std::endl
				<< std::endl;

			//*******************  REMESHING **************************//
			bool remesh = false;
			double ARtol = 1.2;
			uint elementsChanged = 0;

			if (remesh) {
				elementsChanged = bd->Remesh(ARtol, bending, quadOrder);

				//Print out the number of elements that changed due to remeshing
				if (elementsChanged > 0) {
					std::cout << "Number of elements that changed after remeshing = "
						<< elementsChanged << "." << std::endl;

					//We also need to recompute the neighbors for PotentialBody
					PrBody->recomputeNeighbors(PotentialSearchRF);

					//Relax again after remeshing
					solver.solve(&model);

					//Fill-in Matrix for new state for Kabsch algorithm
					for (int col = 0; col < allDefNodes.size(); col++) {
						Vector3D coord = allDefNodes[col]->point();
						for (int row = 0; row < 3; row++) {
							current(row, col) = coord(row);
						}
					}

					rsEnergy = rs.energy();
					bkEnergy = bk.energy();
					vrEnergy = vr.energy();
					bdEnergy = bd->energy() - rsEnergy - vrEnergy - bkEnergy;
					PrEnergy = PrBody->energy();
					energy = solver.function();
				}
			}

			//*********************************************************//

			std::cout << "Shape relaxed." << std::endl
				<< "Energy = " << energy << std::endl;


			//********** Print relaxed configuration ************// 
			Eigen::Affine3d A;
			Eigen::Matrix3Xd newCurr(3, allDefNodes.size());
			A = Find3DAffineTransform(current, initial);
			for (int col = 0; col < current.cols(); col++) {
				newCurr.col(col) = A.linear()*current.col(col)
					+ A.translation();
			}

			//Update the Kabsch transformed positions as new current 
			//configuration in the node container
			for (int i = 0; i < allDefNodes.size(); i++) {
				DeformationNode<3>::Point kabschPoint;
				for (int j = 0; j < 3; j++) {
					kabschPoint(j) = newCurr(j, i);
				}
				allDefNodes[i]->setPoint(kabschPoint);
			}

			bool printQuadPoints = false;

			//We will print only after every currPrintStep iterations
			if (viter % printStep == 0) {
				paraviewStep++;
				sstm << fname << "-relaxed-" << nameSuffix;
				rName = sstm.str();
				//bd->printParaview(rName.c_str());
				bd->printParaview(rName, newCurr, fineConnectivities);
				sstm << ".vtk";
				rName = sstm.str();
				//Store the printed out file name for post-processing
				allVTKFiles.push_back(rName);

				percentStrainData.push_back(percentStrain);

				sstm.str("");
				sstm.clear();
				if (printQuadPoints) {
					sstm << fname << "-QP-" << nameSuffix << ".vtk";
					rName = sstm.str();
					bd->printQuadPoints(rName.c_str());
				}
				sstm.str("");
				sstm.clear();
			}
			//****************************************************//

			//Re-calculate triangle edge length mean and deviation
			lengthStat = calcEdgeLenAndStdDev(allDefNodes, fineConnectivities);
			EdgeLength = lengthStat[0];
			stdDevEdgeLen = lengthStat[1];
			std::cout << "After relaxing with Morse potential: " << endl
				<< "   Average triangle edge length = " << std::setprecision(10)
				<< EdgeLength << endl
				<< "   Standard deviation = " << std::setprecision(10)
				<< stdDevEdgeLen << endl;
			std::cout.precision(6);

			//Calculate centre of sphere as average of position vectors of all nodes.
			tvmet::Vector<double, 3> Xavg(0.0);
			for (int i = 0; i < allDefNodes.size(); i++) {
				Xavg += allDefNodes[i]->point();
			}
			Xavg /= allDefNodes.size();

			//We will calculate radius using the quadrature points
			std::vector<double> radialStats = getRadialStats(bd, Xavg);
			Ravg = radialStats[0];
			std::cout << "Radius of capsid after relaxation = " << Ravg << endl;
			double asphericity = radialStats[1];

			double msd = PrBody->rmsd();
			msd /= (Rshift*Rshift);

			int paraviewStepPrint;
			paraviewStepPrint = (viter % printStep == 0) ? paraviewStep : -1;

			myfile << nameSuffix++ << "\t\t" << paraviewStepPrint << "\t\t"
				<< diffusionCoeff << "\t\t"
				<< Ravg << "\t\t" << asphericity << "\t\t" << gamma << "\t\t"
				<< bdEnergy << "\t\t" << PrEnergy << "\t\t"
				<< rsEnergy << "\t\t" << bkEnergy << "\t\t"
				<< vrEnergy << "\t\t" << energy << "\t\t" << msd
				<< endl;

			//********** Find bins for each particle ************//
			putParticlesInBins(cellLimits, newCurr, defNodes, binDensity, viterMax);

			// step forward in "time", relaxing viscous energy & forces 
			vr.step();
			delete bd;
		}

		pd->GetCellData()->AddArray(binDensity);

		vtkSmartPointer<vtkPolyDataWriter> wr =
			vtkSmartPointer<vtkPolyDataWriter>::New();
		sstm << "density-" << q << ".vtk";
		std::string densityFileName = sstm.str();

		sstm.str("");
		sstm.clear();

		wr->SetFileName(densityFileName.c_str());
		wr->SetInputData(pd);
		wr->Update();
		wr->Write();

		stepCount += viterMax;
	}
	myfile.close();
	t2 = clock();
	float diff((float)t2 - (float)t1);
	std::cout << "Solution loop execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;

	// Post-processing: Manipulating VTK files  
	insertValenceInVtk(allVTKFiles);
	writeEdgeStrainVtk(allVTKFiles, Rshift, percentStrainData);
	t3 = clock();
	diff = ((float)t3 - (float)t2);
	std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
		<< " seconds" << std::endl;
}
