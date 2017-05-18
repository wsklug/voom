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
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"
#include "ViscosityBody.h"

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkGeometryFilter.h>
#include <vtkLinearSubdivisionFilter.h>

#include "Morse.h"
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
  clock_t t1, t2;
  t1 = clock();
  if (argc < 2) {
    cout << "Usage: twoMeshIndent modelName(without .vtk extension)"
	 << endl;
    return(0);
  }

#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  bool verbose = true;
#ifdef WITH_MPI
  MPI_Init(&argc, &argv);
  int procId = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  if (procId != 0) verbose = false;
#endif

  string modelName = argv[1];

  double gamma_inp = 0.0;
  double indent_inp = 0.0;
  double friction_inp = 0.0;
  double step_inp = 0.0;
  bool remesh = false;
  bool unload = false;
  double epsilon;
  double percentStrain;
  double Ravg = 0;
  double Rshift = 1.0;
  double Zmin;
  double Zmax;
  double afmR;
  tvmet::Vector<double, 3> xc(0.0);
  double Z_glass;
  double dZ;
  double ARtol = 1.50;
  double searchRadFactor = 1.2;
  double capsoSearchRadFactor = 1.2;
  double cleanTol = 0.0;
  int numSubDiv = 1;

  //Read epsilon and percentStrain from input file. percentStrain is
  //calculated so as to set the inflection point of Morse potential
  //at a fixed distance relative to the equilibrium separation
  //e.g. 1.1*R_eq, 1.5*R_eq etc.
  string temp;
  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  miscInpFile >> temp >> epsilon
	      >> temp >> percentStrain
	      >> temp >> gamma_inp
	      >> temp >> indent_inp
	      >> temp >> step_inp
	      >> temp >> friction_inp
	      >> temp >> unload
	      >> temp >> ARtol
	      >> temp >> numSubDiv
	      >> temp >> searchRadFactor
	      >> temp >> capsoSearchRadFactor
	      >> temp >> cleanTol;

  miscInpFile.close();

  if (gamma_inp <= 0.0) {
    std::cout << "gamma = " << gamma_inp << " but should be positive."
	      << std::endl;
    return 0;
  }

  vtkPolyDataReader * reader = vtkPolyDataReader::New();
  string inputFileName;

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
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  normals->Update();
  vtkSmartPointer<vtkPolyData> coarseMesh = normals->GetOutput();
  
  //Subdivide into finer mesh
  vtkSmartPointer<vtkLinearSubdivisionFilter> lsdf
    = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  lsdf->SetNumberOfSubdivisions(numSubDiv);
  lsdf->SetInputConnection(normals->GetOutputPort());
  lsdf->Update();
  vtkSmartPointer<vtkPolyData> fineMesh = lsdf->GetOutput();
  
  // create vector of nodes
  int dof = 0;
  std::vector< NodeBase* > bendingNodes;
  std::vector< DeformationNode<3>* > morseNodes;

  std::vector<int> coarseMeshPointIDs(coarseMesh->GetNumberOfPoints(),0);

  // read in points
  for (int a = 0; a < fineMesh->GetNumberOfPoints(); a++) {
    int id = a;
    DeformationNode<3>::Point x, y;
    fineMesh->GetPoint(a, &(x[0]));
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for (int j = 0; j < 3; j++) idx[j] = dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
    bendingNodes.push_back(n);
    for(int b=0; b < coarseMesh->GetNumberOfPoints(); b++){
      coarseMesh->GetPoint(b, &y[0]);
      if(tvmet::norm2(x-y) < 1e-8){
	coarseMeshPointIDs[b] = a;
	morseNodes.push_back(n);
	break;
      }
    }
    
  }

  assert(bendingNodes.size() != 0);
  Ravg /= bendingNodes.size();
  cout << "Number of bending nodes: " << bendingNodes.size() << endl
       << "Ravg = " << Ravg << endl;

  // read in triangle connectivities for fineMesh
  vector< tvmet::Vector<int, 3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri = fineMesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  if (verbose) cout << "Number of triangles in fine mesh: " << ntri << endl;
  for (int i = 0; i < ntri; i++) {
    assert(fineMesh->GetCell(i)->GetNumberOfPoints() == 3);
    for (int a = 0; a < 3; a++) c[a] = fineMesh->GetCell(i)->GetPointId(a);
    connectivities.push_back(c);
  }

  // read in triangle connectivities for coarseMesh
  vector< tvmet::Vector<int, 3> > coarseMeshconn;
  int coarseMeshTri = fineMesh->GetNumberOfCells();
  connectivities.reserve(coarseMeshTri);
  if (verbose) cout << "Number of triangles in coarse mesh: " << coarseMeshTri
		    << endl;
  for (int i = 0; i < coarseMeshTtri; i++) {
    assert(fineMesh->GetCell(i)->GetNumberOfPoints() == 3);
    for (int a = 0; a < 3; a++){
      c[a] = coarseMeshPointIDs[fineMesh->GetCell(i)->GetPointId(a)];
    }
    coarseMeshConn.push_back(c);
  }

  std::vector<double> lengthStat;
  double EdgeLength;
  double stdDevEdgeLen;

  // Calculate side lengths average and std dev of the 
  //equilateral triangles
  lengthStat = calcEdgeLenAndStdDev(bendingNodes, coarseMeshConn);
  EdgeLength = lengthStat[0];
  stdDevEdgeLen = lengthStat[1];
  std::cout << "Before any relaxation :" << endl
	    << "   Average triangle edge length of coarse mesh = "
	    << EdgeLength << endl << "   Standard deviation = "
	    << stdDevEdgeLen << endl;

  //Rescale the capsid such that coarse mesh triangle edge-lengths are
  //unity
  for (int i = 0; i < bendingNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = bendingNodes[i]->point();
    x *= 1.0 / EdgeLength;
    bendingNodes[i]->setPoint(x);
    bendingNodes[i]->setPosition(x);
  }

  //Recalculate edge lengths and capsid radius
  lengthStat = calcEdgeLenAndStdDev(bendingNodes, coarseMeshConn);
  EdgeLength = lengthStat[0];

  Ravg = 0.0;
  for (int i = 0; i < bendingNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = bendingNodes[i]->point();
    double tempRadius = tvmet::norm2(x);
    Ravg += tempRadius;
  }
  Ravg /= bendingNodes.size();

  std::cout << "Radius of capsid after rescaling = " << Ravg << endl;

  //Set the AFM indenter radius
  afmR = Ravg;

  //Material properties
  Rshift = EdgeLength;
  double sigma = (100 / (Rshift*percentStrain))*log(2.0);
  double PotentialSearchRF = searchRadFactor*Rshift;
  double springConstant = 2 * sigma*sigma*epsilon;
  double gamma = gamma_inp;
  double Y = 2.0 / sqrt(3)*springConstant; //2D Young's modulus
  double nu = 1.0 / 3.0;
  double KC = Y*Ravg*Ravg / gamma; //Bending modulus
  double KG = -2 * (1 - nu)*KC; // Gaussian modulus
  double C0 = 0.0;
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
  MaterialType bending(KC, KG, C0, 0.0, 0.0);

  //********************************** Actual Indentation **********************//

  LSB * bd = new LSB(bending, connectivities, bendingNodes, quadOrder);
  bd->setOutput(paraview);

  std::cout << "Morse  potential parameters:" << endl
	    << "\tsigma = " << sigma << endl
	    << "\tPotentialSearchRF = " << PotentialSearchRF << endl
	    << "\tepsilon = " << epsilon << endl
	    << "\tIs Remesh On? " << remesh << endl
	    << "\tARtol = " << ARtol << endl;
		
  // Protein body implemented using Morse potential body
  Morse Mat(epsilon, sigma, Rshift);

  // Then initialize potential body
  PotentialBody * PrBody = new PotentialBody(&Mat, morseNodes, PotentialSearchRF);
  PrBody->compute(true, false, false);
  std::cout << "Initial protein body energy = " << PrBody->energy() << endl;
  //Assign the MorseBond structure from the initial neighbor information
  vtkSmartPointer<vtkCellArray> bonds = vtkSmartPointer<vtkCellArray>::New();
  getMorseBonds(bonds, morseNodes, PotentialSearchRF);

  // Create a ViscosityBody for viscous regularization
  ViscosityBody * vb = new ViscosityBody(fineMesh, connectivities,
					 1e-2*springConstant);
  
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(PrBody);
  bdc.push_back(bd);
  bdc.push_back(vb);

  Model model(bdc, nodes);

  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, 1e5);
  std::cout << "Relaxing shape for gamma = " << gamma << std::endl;

  // relax initial shape;
  while(vb->energy() > 1e-8*(bd->energy() + PrBody->energy())){
    solver.solve(&model);
  }

  std::cout << "Shape relaxed." << std::endl
	    << "Energy = " << solver.function() << std::endl;

  //Calculate centre of sphere as average of position vectors of all nodes.
  tvmet::Vector<double, 3> Xavg(0.0);
  for (int i = 0; i < bendingNodes.size(); i++) {
    Xavg += bendingNodes[i]->point();
  }
  Xavg /= bendingNodes.size();

  //We will calculate radius using the quadrature points
  vtkSmartPointer<vtkPolyData> lssPd = bd->getLoopShellSurfPoints(cleanTol);
  std::vector<double> radialStats = getRadialStats(lssPd, Xavg);
  Ravg = radialStats[0];
  std::cout << "Radius of capsid after relaxation = " << Ravg << endl;
  double asphericity = radialStats[1];
  double capsomerSearchRad = radialStats[2];

  std::cout << "Effective 2D Young's modulus = " << Y << endl
	    << "FVK number = " << gamma << endl
	    << "Asphericity = " << asphericity << endl;

  // find top and bottom of capsid
  Zmin = std::numeric_limits<double>::max();
  Zmax = -std::numeric_limits<double>::max();

  for (int a = 0; a < bendingNodes.size(); a++) {
    double Z = bendingNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
    Zmax = std::max(Zmax, Z);
  }

  std::cout << "Zmax = " << Zmax << std::endl
	    << "Zmin = " << Zmin << std::endl;

  xc = 0.0, 0.0, Zmax + afmR;

  std::cout << "AFM Indenter radius =" << afmR << std::endl
	    << "AFM Indenter center = (" << xc[0] << ","
	    << xc[1] << "," << xc[2] << ")" << std::endl;

  std::cout << "Compressing capsid..........." << std::endl;

  double friction = friction_inp;
  double k_AL = 1.0e2;
  RigidHemisphereAL * afm
    = new RigidHemisphereAL(bendingNodes, k_AL, afmR, xc, friction);
  afm->updateContact();

  bd->pushBack(afm);
  std::cout << "Added afm to body." << std::endl;

  bool up = true;
  Z_glass = Zmin;

  RigidPlateAL* glass = new RigidPlateAL(bendingNodes, k_AL, Z_glass, up, friction);
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
       << "\tMorseEnergy\tLSBStrainEnergy\tPlateEnergy\tAFMEnergy"
       <<"\tViscRegEnergy\tSolverFunction"
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

  // Following variables is for output from LoopShellBody::Remesh()
  uint elementsChanged = 0;

  for (double Z = Zbegin; ; Z += dZ, step++) {

    // initial guess
    if (step == 0) {
      // shift capsid up by dZ/2 as an initial guess
      for (int a = 0; a < bendingNodes.size(); a++) {
	bendingNodes[a]->addPoint(2, 0.5*dZ);
      }
      model.getField(solver);
    }
    else if (std::abs(dZ) > 0) { 
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
      PrBody->checkConsistency(true,false);
      }
    */

    model.computeAndAssemble(solver, false, true, false);
    bd->printParaview("contact");

    while(vb->energy() > 1e-8*(bd->energy() + PrBody->energy())){
      solver.solve(&model);
    }
    
    bdEnergy = bd->energy();		

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
    double Ztop = Zmax;
    double Zbot = Zmax - height;

    if (verbose) {
      std::cout << "Contact converged." << std::endl
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
    if (max(abs(u_prev)) > 1.0) {
      std::cout << "??????? Large Deformation Detected ???????"
		<< std::endl;
      u_prev = 0.0;//????
    }

    // Save current (successful) state as previous
    for (int i = 0; i < model.dof(); i++) x_prev(i) = solver.field(i);

    // reached max indentation, now unload by reversing dZ
    if (unload && Z >= Zend) { 
      dZ = -dZ;
    }

    //This is where we remesh before next indentation step
    if (remesh) {
      elementsChanged = bd->Remesh(ARtol, bending, quadOrder);
			
      //Print out the number of elements that changed due to remeshing
      if (elementsChanged > 0) {
	std::cout << "Number of elements that changed after remeshing = "
		  << elementsChanged << "." << std::endl;
	
	//Relax again after remeshing
	while(vb->energy() > 1e-8*(bd->energy() + PrBody->energy())){
	  solver.solve(&model);
	}

      }
    }// Remeshing ends here

    //*********** BEGIN PRINTING OUTPUT (and log) FILES ***********//

    FvsZ << step
	 << "\t" << (originalHeight - height) / Ravg
	 << "\t" << glass->FZ() / std::sqrt(Y*KC)
	 << "\t" << afm->FZ() / std::sqrt(Y*KC)
	 << "\t" << PrBody->energy() / (Ravg*std::sqrt(Y*KC))
	 << "\t" << bd->totalStrainEnergy() / (Ravg*std::sqrt(Y*KC))
	 << "\t" << glass->energy() / (Ravg*std::sqrt(Y*KC))
	 << "\t" << afm->energy() / (Ravg*std::sqrt(Y*KC))
	 << "\t" << vb->energy() / (Ravg*std::sqrt(Y*KC))
	 << "\t" << solver.function() / (Ravg*std::sqrt(Y*KC))
	 << std::endl;

    sstm << modelName << "-step-" << step;
    rName = sstm.str();
    bd->printParaview(rName);
    //We will append Caspsomer POINT_DATA to the vtk output file
    //printed by printParaview(), if such a file exists
    sstm << ".vtk";
    rName = sstm.str();
    sstm.str("");
    sstm.clear();		
    allStepFiles.push_back(rName);
			
    //Now we will print the LoopShellSurface
    //We will calculate radius using the quadrature points
    lssPd = bd->getLoopShellSurfPoints(cleanTol);
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
  //writeEdgeStrainVtk(fakeVec, Rshift, percentStrain);
  plotMorseBonds(allStepFiles, modelName, epsilon, Rshift, sigma, bonds);
  FvsZ.close();

  std::cout << "Indentation complete." << std::endl;

  t2 = clock();
  float diff((float)t2 - (float)t1);
  std::cout << "Total execution time: " << diff / CLOCKS_PER_SEC
	    << " seconds" << std::endl;
  return 0;

}
