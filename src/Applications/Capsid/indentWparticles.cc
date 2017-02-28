#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"

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
  clock_t t1, t2;
  t1 = clock();
  if (argc < 2) {
    cout << "Usage: indent modelName"
	 << endl;
    return(0);
  }

#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  bool verbose = false;
#ifdef WITH_MPI
  MPI_Init(&argc, &argv);
  int procId = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);
  if (procId != 0) verbose = false;
#endif

  for (int i = 0; i < argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

  string modelName = argv[1];

  double gamma_inp = 0.0;
  double indent_inp = 0.0;
  double viscosity_inp = 0.0;
  double friction_inp = 0.0;
  double step_inp = 0.0;
  bool remesh = true;
  bool unload = false;

  std::ifstream cmdInp("cmdInp.dat");
  assert(cmdInp);
  string temp;
  cmdInp >> temp >> gamma_inp
	 >> temp >> indent_inp
	 >> temp >> step_inp
	 >> temp >> friction_inp
	 >> temp >> viscosity_inp
	 >> temp >> unload;

  if (gamma_inp <= 0.0) {
    std::cout << "gamma = " << gamma_inp << " but should be positive."
	      << std::endl;
    return 0;
  }

  //For Morse material 
  double epsilon;
  double percentStrain;
  double pressureFactor;
  bool harmonicRelaxNeeded;
  double Ravg = 0;
  double Rshift = 1.0;
  double Zmin;
  double Zmax;
  double afmR;
  tvmet::Vector<double, 3> xc(0.0);
  double Z_glass;
  double dZ;

  //Read epsilon and percentStrain from input file. percentStrain is
  //calculated so as to set the inflection point of Morse potential
  //at a fixed distance relative to the equilibrium separation
  //e.g. 1.1*R_eq, 1.5*R_eq etc.
  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  miscInpFile >> temp >> epsilon
	      >> temp >> percentStrain
	      >> temp >> pressureFactor
	      >> temp >> harmonicRelaxNeeded;

  miscInpFile.close();

  vtkPolyDataReader * reader = vtkPolyDataReader::New();
  vtkSmartPointer<vtkPolyData> mesh;
  string inputFileName;
  vtkSmartPointer<vtkDataArray> displacements;

  vtkSmartPointer<vtkPolyData> mesh_prev;
  vtkSmartPointer<vtkDataArray> displacements_prev;

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
  Rshift = EdgeLength;
  double sigma = (100 / (Rshift*percentStrain))*log(2.0);
  double PotentialSearchRF = 1.2*Rshift;
  double springConstant = 2 * sigma*sigma*epsilon;
  double pressure = 12 * sigma*epsilon
    *(exp(-2 * sigma*Rshift) - exp(-sigma*Rshift)
      + exp(-1.46410*sigma*Rshift) - exp(-0.7321*sigma*Rshift))
    / (3 * Ravg*Ravg);

  if (pressure < 0.0) {
    pressure = pressure*(-1);
  }

  double fracturePressure = (3.82)*sigma*epsilon / (Rshift*Rshift);
  std::cout << "Fracture Pressure = " << fracturePressure << endl
	    << "Minimum Pressure = " << pressure << endl;
  pressure *= pressureFactor;
  std::cout << "Pressure in use = " << pressure << endl;

  double gamma = gamma_inp;
  double Y = 2.0 / sqrt(3)*springConstant; //2D Young's modulus
  double nu = 1.0 / 3.0;
  double KC = Y*Ravg*Ravg / gamma; //Bending modulus
  double KG = -2 * (1 - nu)*KC; // Gaussian modulus
  C0 = 0.0;
  double ARtol = 1.1;
  int quadOrder = 2;

  int m = 5;
  int maxIter = 1e5;
  double factr = 1.0e+1;
  double pgtol = 1.0e-7;
  int iprint = 1;

  std::stringstream sstm;
  string fname = modelName;
  string rName;
  string actualFile;

  typedef FVK MaterialType;
  typedef LoopShellBody<MaterialType> LSB;
  typedef LoopShell<MaterialType> LS;
  MaterialType bending(KC, KG, C0, 0.0, 0.0);

  if (harmonicRelaxNeeded) {

    //****** Relax the initial mesh using harmonic potential ****** //

    LSB * bd1 = new LSB(bending, connectivities, nodes, quadOrder,
			pressure, 0.0, 0.0, 1.0e4, 1.0e6, 1.0e4,
			multiplier, noConstraint, noConstraint);
    bd1->setOutput(paraview);
    SpringPotential SpringMat(springConstant, Rshift);
    PotentialBody * SpringBody = new
      PotentialBody(&SpringMat, defNodes, PotentialSearchRF);

    //Create Model
    Model::BodyContainer bdc1;
    bdc1.push_back(SpringBody);
    bdc1.push_back(bd1);
    Model model1(bdc1, nodes);

    std::cout << "Spring constant: " << springConstant << endl;
    std::cout << "Relaxing the mesh using harmonic potential..." << endl;

    Lbfgsb solver2(model1.dof(), m, factr, pgtol, iprint, 1e5);
    solver2.solve(&model1);

    std::cout << "Harmonic potential relaxation complete." << endl;

    //Print to VTK file
    sstm << fname << "-relaxed-harmonic";
    rName = sstm.str();
    bd1->printParaview(rName.c_str());
    sstm.str("");
    sstm.clear();
    
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
	      << "   Average triangle edge length = " << std::setprecision(10)
	      << EdgeLength << endl
	      << "   Standard deviation = " << std::setprecision(10)
	      << stdDevEdgeLen << endl;
    std::cout.precision(6);
  }

  //********************************** Actual Indentation **********************//

  LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
  bd->setOutput(paraview);

  std::cout << "Morse  potential parameters:" << endl
	    << "sigma = " << sigma << endl
	    << "PotentialSearchRF = " << PotentialSearchRF << endl
	    << "epsilon = " << epsilon << endl
	    << "Is Remesh On? " << remesh << endl
	    << "ARtol = " << ARtol << endl;

  // Protein body implemented using Morse potential body
  Morse Mat(epsilon, sigma, Rshift);

  // Then initialize potential body
  PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);
  PrBody->compute(true, false, false);
  std::cout << "Initial protein body energy = " << PrBody->energy() << endl;

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(PrBody);
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
  std::vector<double> radialStats = getRadialStats(bd, Xavg);
  Ravg = radialStats[0];
  std::cout << "Radius of capsid after relaxation = " << Ravg << endl;
  double asphericity = radialStats[1];

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
 
  // add some viscosity for regularization
  double minViscosity = 1.0e-6*viscosity_inp;
  double maxViscosity = 1.0e+6*viscosity_inp;
  ViscousRegularizer vr(bd->nodes(), viscosity_inp);
  bd->pushBack(&vr);

  // set viscosity parameters
  double targetVelocity = std::abs(dZ);
  double vrTol = 1.0e-10;

  double vrEnergy = vr.energy();
  double bdEnergy = bd->energy();

  string fzName = modelName + ".fz";
  ofstream FvsZ(fzName.c_str());
  FvsZ << "#Step\tIndentation\tGlass_Fz\tAFM_Fz"
       <<"\tElasticEnergy\tMorseEnergy"
       <<std::endl;

  blitz::Array<double, 1> x_prev(model.dof());
  blitz::Array<double, 1> u_prev(model.dof());
 
  model.getField(solver);
  for (int i = 0; i < model.dof(); i++) x_prev(i) = solver.field(i);
  u_prev = 0.0;
  
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

    int viterMax = 20;
    for (int viter = 0; viter < viterMax; viter++) {
      if (viter == viterMax - 1) vr.setViscosity(minViscosity);
      blitz::Array<double, 1> vSave(model.dof());
      model.getField(solver);
      for (int i = 0; i < model.dof(); i++) vSave(i) = solver.field(i);

      if (verbose) std::cout << std::endl
			     << "VISCOUS ITERATION: " << viter
			     << "\t viscosity = " << vr.viscosity()
			     << std::endl
			     << std::endl;
      // update contact
      glass->updateContact();
      afm->updateContact();

      model.computeAndAssemble(solver, false, true, false);
      bd->printParaview("contact");

      solver.solve(&model);
      vrEnergy = vr.energy();
      bdEnergy = bd->energy();

      if (verbose) {
	std::cout << "ENERGY:" << std::endl
		  << "viscous energy = " << vrEnergy << std::endl
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

      if (verbose) {
	std::cout << "VISCOSITY: " << std::endl
		  << "          velocity = " << vr.velocity() << std::endl
		  << "   target velocity = " << targetVelocity << std::endl
		  << " updated viscosity = " << vr.viscosity() << std::endl
		  << std::endl;
      }

      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if (vrEnergy < std::abs(vrTol*bdEnergy) && solver.projectedGradientNorm() <= pgtol &&
	  std::max(afm->penetration(), glass->penetration()) < 1.0e-2*std::abs(dZ)) {
	// viscous energy is small enough; exit
	break;
      }
    }

    double height = Zmax - Z;
    // add up forces on top and bottom
    double Ztop = Zmax;
    double Zbot = Zmax - height;

    if (verbose) {
      std::cout << "Contact and Viscous energy converged." << std::endl
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

    //This is where we remesh before next indentation step

    if (remesh) {
      elementsChanged = bd->Remesh(ARtol, bending, quadOrder);

      //Print out the number of elements that changed due to remeshing
      if (elementsChanged > 0) {
	std::cout << "Number of elements that changed after remeshing = "
		  << elementsChanged << "." << std::endl;

	//If some elements have changed then we need to reset the
	//reference configuration with average side lengths
	bd->SetRefConfiguration(EdgeLength);

	//We also need to recompute the neighbors for PotentialBody
	PrBody->recomputeNeighbors(PotentialSearchRF);

	//Relax again after remeshing
	solver.solve(&model);

      }
    }// Remeshing ends here

    //*********** BEGIN PRINTING OUTPUT (and log) FILES ***********//

    FvsZ << std::setw(10)
	 << step
	 << std::setw(24) << std::setprecision(16)
	 << (originalHeight - height)/Ravg
	 << std::setw(24) << std::setprecision(16)
	 << glass->FZ()/std::sqrt(Y*KC)
	 << std::setw(24) << std::setprecision(16)
	 << afm->FZ()/std::sqrt(Y*KC)
	 << std::setw(24) << std::setprecision(16)
	 << bd->energy()/(Ravg*std::sqrt(Y*KC))
	 << std::setw(24) << std::setprecision(16)
	 << PrBody->energy()/(Ravg*std::sqrt(Y*KC))
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
    if (ifstream(rName.c_str())) {
      std::vector<string> fakeVec;
      fakeVec.push_back(rName);
      insertValenceInVtk(fakeVec);
      writeEdgeStrainVtk(fakeVec, Rshift, percentStrain);
    }
    //************* END PRINTING OUTPUT FILES **************//

    // check if we are done
    if (unload && Z + dZ < Zbegin)
      break;
    else if (!unload && Z + dZ > Zend)
      break;

  }// Indentation Loop Ends

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
