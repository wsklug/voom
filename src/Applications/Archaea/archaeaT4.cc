#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

#include <getopt.h>
#include <unistd.h>
#include <time.h>

#include <tvmet/Vector.h>

#include "Node.h"
#include "EvansElastic.h"
#include "LoopShellBody.h"
#include "PotentialBody.h"
// #include "SpringPotential.h"
// #include "SpringPotentialSQ.h"
// #include "LennardJones.h"
// #include "LennardJonesFT.h"
#include "Morse.h"
#include "ViscosityBody.h"
#include "LoopShell.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "VoomMath.h"
#include "Utils/PrintingArchaeaS.h"

using namespace voom;
using namespace std;

int main(int argc, char* argv[])
{

  // -------------------------------------------------------------------
  // Setup
  // -------------------------------------------------------------------
  time_t start, end;
  double dif;
  time(&start);

  if(argc < 2){
    cout << "Input file missing." << endl;
    return(0);
  }

  // File names
  string parameterFileName = argv[1];
  string modelName;
  string outputFileName;
  uint restart = 0;
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;

  // Applied force parameters
  double ForceSearchR = 0.0;
  double MaxDisp = 1.0;
  double nDispSteps = 1.0;

  // Regularization algorithm parameters
  double kvisc = 1.0;
  double viscTol = 1.0e-4;

  // Penalty parameters
  double penaltyVolumeInit = 1.0e3;
  double penaltyAreaInit   = 1.0e3;
  double penaltyVolumeTol  = 1.0e-3;
  double penaltyAreaTol    = 1.0e-3;

  // Potential input parameters
  double PotentialSearchRF = 1.0;
  double epsilon = 1.0;
  double sigma = 1.0;
  double Rshift = 1.0;

  // Connectivity search
  double RconnSF = 1.0;

  // Remeshing parameters
  double ARtol = 2.0;
  
  // Reduced volume
  double Vmin = 0.9;

  // PotentialBodyFull flag
  uint PotentialBodyFull = 1;  // 0 -> nodes at which we apply displacements are not proteins (no Morse potential interactions)
                               // 1 -> nodes at which we apply displacements are proteins (Included in potential body)

  int quadOrder = 1;

  // Reading input from file passed as argument
  ifstream inp;
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  string temp;
  
  inp >> temp >> modelName;
  inp >> temp >> outputFileName;
  inp >> temp >> restart;
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> ForceSearchR;
  inp >> temp >> MaxDisp;
  inp >> temp >> nDispSteps;
  inp >> temp >> kvisc;
  inp >> temp >> viscTol;
  inp >> temp >> penaltyVolumeInit;
  inp >> temp >> penaltyAreaInit;
  inp >> temp >> penaltyVolumeTol;
  inp >> temp >> penaltyAreaTol;
  inp >> temp >> PotentialSearchRF;
  inp >> temp >> epsilon;
  inp >> temp >> sigma;
  inp >> temp >> Rshift;
  inp >> temp >> RconnSF;
  inp >> temp >> ARtol;
  inp >> temp >> Vmin;
  inp >> temp >> PotentialBodyFull;
  inp >> temp >> quadOrder;


  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " outputFileName          : " << outputFileName    << endl
       << " restart                 : " << restart           << endl
       << " KC                      : " << KC                << endl
       << " KG                      : " << KG                << endl
       << " C0                      : " << C0                << endl
       << " Force Search Radius     : " << ForceSearchR      << endl
       << " Max Disp                : " << MaxDisp           << endl
       << " N disp steps            : " << nDispSteps        << endl
       << " Viscous regularization  : " << kvisc             << endl
       << " Viscous reg tolerance   : " << viscTol           << endl
       << " Penalty volume factor   : " << penaltyVolumeInit << endl
       << " Penalty area factor     : " << penaltyAreaInit   << endl  
       << " Penalty volume tol      : " << penaltyVolumeTol  << endl
       << " Penalty area tol        : " << penaltyAreaTol    << endl  
       << " Potential Search factor : " << PotentialSearchRF << endl
       << " epsilon                 : " << epsilon           << endl
       << " sigma                   : " << sigma             << endl
       << " Rshift                  : " << Rshift            << endl
       << " Rconn scaling factr     : " << RconnSF           << endl
       << " Aspect Ratio tol        : " << ARtol             << endl
       << " Min reduced volume      : " << Vmin              << endl
       << " PotentialBodyFull       : " << PotentialBodyFull << endl
       << " Loop Shell quadOrder    : " << quadOrder         << endl;


      
  






  // -------------------------------------------------------------------
  // Read mesh and setup geometric model
  // -------------------------------------------------------------------
  // Read in the mesh, in (ascii) legacy vtk format.
  // Create input stream
  ifstream ifs;
  ifs.open(modelName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << modelName << endl;
    exit(0);
  }
 
  // Create vector of nodes
  unsigned int dof = 0, npts = 0, NumDoF = 0;
  vector<NodeBase* > nodes;
  vector<DeformationNode<3>* > defNodes, defNodesV0, defNodesV1, potentialNodes;
  double Ravg = 0;

  // Opposite vertices where forces are applied
  Vector3D V0(0.0, 0.0, 1.902113);
  Vector3D V1(0.0, 0.0,-1.902113);
  uint CountV0nodes = 0, CountV1nodes = 0, CountSolverNodes = 0;

  // Bounds for solver
    blitz::Array<double,1> LowerBound, UpperBound;
    blitz::Array<int,1> BoundType;

   
  
  // Input .vtk or inp files containing nodes and connectivities
  string token;

  ifs >> token;
  while( token != "POINTS" ) ifs >> token;
  ifs >> npts; 
  NumDoF = npts*3; // Assumed all nodes have 3 DoF
  nodes.reserve(npts);
  defNodes.reserve(npts);

  // Initialize bound of dof
  LowerBound.resize(NumDoF);
  UpperBound.resize(NumDoF);
  BoundType.resize(NumDoF);
  LowerBound = 0.0;
  UpperBound = 0.0;
  BoundType  = 0;

  ifs >> token;   // skip number type in vtk file
  // read in points
  for(uint i = 0; i < npts; i++) {
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) idx[j] = dof++;
    
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
    
    // Distinguish between free and prescribed nodes
    // Check only distance in x,y plane so can use it also when restart = true
    if ( (sqrt(pow(x(0)-V0(0), 2.0) + pow(x(1)-V0(1), 2.0)) < ForceSearchR) && (x(2) > 0.0) ) {
      defNodesV0.push_back(n);
    }
    else if ( (sqrt(pow(x(0)-V1(0), 2.0) + pow(x(1)-V1(1), 2.0)) < ForceSearchR) && (x(2) < 0.0) )  {
      defNodesV1.push_back(n);
    }
    else if (PotentialBodyFull == 0) {
      potentialNodes.push_back( n );
    }
   
  }
  Ravg /= nodes.size();
  
  CountV0nodes = defNodesV0.size();
  CountV1nodes = defNodesV1.size();
  
  cout << endl << "V0 nodes number = " << CountV0nodes << endl;
  cout << endl << "V1 nodes number = " << CountV1nodes << endl;
  // assert(CountV0nodes == CountV1nodes);
  

   
  // Read in triangle connectivities
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  uint ntri = 0, ElemNum = 0, tmp = 0;
  double AverageEdgeLength = 0.0;

  while( token != "POLYGONS" ) ifs >> token;
  ifs >> ntri;
  connectivities.reserve(ntri);
  ifs >> tmp;
  for (uint i = 0; i < ntri; i++)
  {
    ifs >> tmp;
    if(tmp != 3) cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    ifs >> ct(0) >> ct(1) >> ct(2);
    // Compute initial average edge length - used in potential body
    AverageEdgeLength += tvmet::norm2(defNodes[ct(0)]->point()- defNodes[ct(1)]->point());
    AverageEdgeLength += tvmet::norm2(defNodes[ct(1)]->point()- defNodes[ct(2)]->point());
    AverageEdgeLength += tvmet::norm2(defNodes[ct(2)]->point()- defNodes[ct(0)]->point());
    
    connectivities.push_back(ct);
  }
 
  AverageEdgeLength /= double(ntri*3);
  
  cout << "Number of nodes     = " << nodes.size()          << endl
       << "Number of elements  = " << connectivities.size() << endl
       << "Ravg                = " << Ravg                  << endl
       << "Average Edge Length = " << AverageEdgeLength     << endl;

  // Close mesh file in vtk or inp format
  ifs.close();
  


  // -------------------------------------------------------------------
  // Create model
  // -------------------------------------------------------------------
  // Create Loop shell body with both bending and in-plane elastic energy
  
  // Bending material object, to be copied when body generates new elements
  EvansElastic bodyMaterial(KC, KG, C0, 0.0, 0.0);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  double pressure = 0.0;
  double tension  = 0.0;
  LoopShellBody<EvansElastic> * shellBody = new LoopShellBody<EvansElastic>(bodyMaterial, 
									    connectivities, 
									    nodes, 
									    quadOrder,
									    pressure, 
									    tension,
									    0.0,
									    penaltyVolumeInit,
									    penaltyAreaInit,
									    0.0,
									    noConstraint, // just at the beginning
									    noConstraint, // just at the beginning
									    noConstraint);
  // Initialize body
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  cout << "Initial shell body volume = " << shellBody->volume() << endl;
  cout << "Initial shell body area   = " << shellBody->area()   << endl;
  cout << "Average edge length from body = " << shellBody->AverageEdgeLength() << endl << endl;
  // shellBody->checkConsistency();
  


  // Create viscosity body
  ViscosityBody RegularizationBody(defNodes, connectivities, 0.0); // epsilon = 0 at the beginning
  RegularizationBody.compute(true, false, false);
  cout << "Initial regularization body energy = " << RegularizationBody.totalStrainEnergy() << endl;
  cout << "Initial regularization body energy = " << RegularizationBody.energy() << endl;



  // Initiliaze potential body
	// First initialize material
	   // double KspringScaled = epsilon/pow(AverageEdgeLength, 4.0);
           // SpringPotentialSQ MatSpring(KspringScaled, LrestScaled);
           // double LrestScaled = sigma*AverageEdgeLength;
	double RshiftScaled = Rshift*AverageEdgeLength;
	Morse ProteinPotential(epsilon, sigma, RshiftScaled); // Forces from protein potential also during deflating

	// Then initialize potential body
  	double ScaledSearchR = PotentialSearchRF*AverageEdgeLength;
	PotentialBody * ProteinBody;
	if (PotentialBodyFull == 1) {
	  ProteinBody = new PotentialBody(&ProteinPotential, defNodes, ScaledSearchR); }
	else {
	  ProteinBody = new PotentialBody(&ProteinPotential, potentialNodes, ScaledSearchR);
	}
	  
  	// potentialBody->checkConsistency();
  
  	ProteinBody->compute(true, false, false);
  	cout << "Initial protein body (FirstRing) energy = " << ProteinBody->totalStrainEnergy() << endl;
  	cout << "Initial protein body (FirstRing) energy = " << ProteinBody->energy() << endl;







  // Create Model and solver
  Model::BodyContainer bdc;
  bdc.push_back(shellBody);
  double TolToAddProteinBody = 1.0e-12;
  if (epsilon > TolToAddProteinBody) {
    bdc.push_back(ProteinBody);
  }
  else {
    cout << "No protein body added." << endl;
  }
  bdc.push_back(&RegularizationBody);
  
  Model model(bdc, nodes);
       // Consistency check
       // model.checkConsistency(true,false);
       // model.checkRank(model.dof()-3,true);
  
  // Set solver and solve
  int m = 5;  
  double factr = 1.0e1;
  double pgtol = 1.0e-8;
  int iprint = 0; 
  int maxIter = 30000;
  cout << endl << "Input iprint: " << iprint << endl
               << "Input factr:  " << factr  << endl
               << "Input pgtol:  " << pgtol  << endl
               << "Input m:      " << m      << endl;

  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);
  cout << "Solver DOF = " << model.dof() << endl;







  // Initialize printing utils
  cout << "Seraching connectivity over R = " <<  double(RconnSF*AverageEdgeLength) << endl;
  PrintingArchaeaS PrintVirus(modelName, outputFileName+"iter", bdc, double(RconnSF*AverageEdgeLength), defNodes, &ProteinPotential);

  // Print before everything
  uint PrintCount = 0;
  vector<uint > Dislocation(3,0);
  Dislocation = PrintVirus.printMaster(PrintCount);
  char vtkname[100]; 
  sprintf(vtkname,"Capsid-%04d",PrintCount);
  shellBody->printParaview(vtkname);

  

  /*
  // Check on remeshing
  uint ElementsChanged = shellBody->Remesh(ARtol, bodyMaterial, quadOrder);
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  cout << "Initial shell body volume = " << shellBody->volume() << endl;
  cout << "Initial shell body area   = " << shellBody->area()   << endl;
  cout << "Average edge length from body = " << shellBody->AverageEdgeLength() << endl;
  cout << "Elements changed = " <<  ElementsChanged << endl << endl;

  PrintVirus.printMaster(++PrintCount);
  sprintf(vtkname,"Capsid-%04d",PrintCount);
  shellBody->printParaview(vtkname);
  */







  // -------------------------------------------------------------------
  // Solve
  // -------------------------------------------------------------------

  // Add constraint to keep deformation symmetric during deflation
  double tol = 1.0e-5;
  for (uint i = 0; i < defNodes.size(); i++)
  {
    DeformationNode<3>::Point x;
    NodeBase::DofIndexMap idx(3);
    x = defNodes[i]->point();
    // if ( fabs(x[0]) > (1.618034-tol) && fabs(x[1]) < tol && fabs(x[2]) < tol) {
    //   idx = defNodes[i]->index();

    //   BoundType(idx[1]) = 2;
    //   LowerBound(idx[1]) = 0.0;
    //   UpperBound(idx[1]) = 0.0;

    //   BoundType(idx[2]) = 2;
    //   LowerBound(idx[2]) = 0.0;
    //   UpperBound(idx[2]) = 0.0;
    // }
    // else if ( fabs(x[0]) < tol && fabs(x[1]) > (1.538842-tol) && fabs(x[2]) < tol) {
    //   idx = defNodes[i]->index();
    //   BoundType(idx[0]) = 2;
    //   LowerBound(idx[0]) = 0.0;
    //   UpperBound(idx[0]) = 0.0;
  
    //   BoundType(idx[2]) = 2;
    //   LowerBound(idx[2]) = 0.0;
    //   UpperBound(idx[2]) = 0.0;
    // }
    // else if ( fabs(x[0]) < tol && fabs(x[1]) < tol && fabs(x[2]) > (1.902113-tol) ) {
    if ( fabs(x[0]) < tol && fabs(x[1]) < tol ) {
      idx = defNodes[i]->index();
      BoundType(idx[0]) = 2;
      LowerBound(idx[0]) = 0.0;
      UpperBound(idx[0]) = 0.0;
  
      BoundType(idx[1]) = 2;
      LowerBound(idx[1]) = 0.0;
      UpperBound(idx[1]) = 0.0;
      // if ( x[2] <  0.0 )
      // {
      // 	BoundType(idx[2]) = 2;
      // 	LowerBound(idx[2]) = x[2];
      // 	UpperBound(idx[2]) = x[2];
      // }
    }
  }



  // Solve to minimize the energy - Initial relaxed configuration
  double FixedShellArea = 0.0;
  double FixedShellVolume = 0.0;
  if (epsilon > TolToAddProteinBody) {
    solver.setBounds(BoundType, LowerBound, UpperBound);
    solver.solve(&model);
    Dislocation = PrintVirus.printMaster(++PrintCount);
    sprintf(vtkname,"Capsid-%04d",PrintCount);
    shellBody->printParaview(vtkname);
    FixedShellArea = shellBody->area();
    FixedShellVolume = shellBody->volume();

       // Compute Force vs Displacement
       Vector3D Force0(0.0), Force1(0.0);
       for (uint i = 0; i < defNodesV0.size(); i++) {
	 Force0 += defNodesV0[i]->force();
       }
       for (uint i = 0; i < defNodesV1.size(); i++) {
	 Force1 += defNodesV1[i]->force();
       }
       cout << "     u = " << defNodesV0[0]->getPoint(2) << "   "
	    << "Force0 = " << Force0(0) << " " << Force0(1) << " " <<Force0(2) << "   "
	    << "Force1 = " << Force1(0) << " " << Force1(1) << " " <<Force1(2) << endl;
       cout << "Pressure = " << shellBody->pressure() << " Surface Tension = " << shellBody->tension() << endl;
       cout << "Vertices valence [5, 6, 7] = " << Dislocation[0] << " " << Dislocation[1] << " " << Dislocation[2] << endl;

  }
  else {
    FixedShellArea = 33.0057;
    FixedShellVolume = 49.9053;
  }

  cout << "Shell body volume before reduction = " << shellBody->volume() << endl;
  cout << "Shell body area before reduction   = " << shellBody->area()   << endl;

  // Now set area to the relaxed value computed using bending energy and protein potential 
  shellBody->setAreaConstraint(augmented);
  shellBody->setPrescribedArea(FixedShellArea);

  if (Vmin < 1.0) {
    shellBody->setVolumeConstraint(augmented);
    shellBody->setPrescribedVolume(FixedShellVolume);
  }

  // Activate Regularization Body too
  RegularizationBody.setSigma(kvisc);
  RegularizationBody.resetRefConf();
 
  // Initiliaze shell by reducing its volume to 0.9
  double Vred = 1.0;


  
  // Now change bounds on nodes on z axis
  for (uint i = 0; i < defNodes.size(); i++)
  {
    DeformationNode<3>::Point x;
    NodeBase::DofIndexMap idx(3);
    x = defNodes[i]->point();
    if ( fabs(x[0]) < tol && fabs(x[1]) < tol ) {
      idx = defNodes[i]->index();
      if ( x[2] > tol ) {
  	BoundType(idx[2]) = 1;
  	LowerBound(idx[2]) = x[2];
      }
      else if ( x[2] < tol) {
  	BoundType(idx[2]) = 3;
  	UpperBound(idx[2]) = x[2];
      }
    }
  }



  if (!restart)
  {
    for(Vred = 0.99; Vred > (Vmin-1.0e-5); Vred -= 0.01)
    {
      shellBody->reduceVolume(Vred);
     
      // Update neighbors used in ProteinBody
      if (epsilon > TolToAddProteinBody) {
	ProteinBody->recomputeNeighbors(ScaledSearchR);
      }
  
      double penaltyVolume = penaltyVolumeInit;
      double penaltyArea   = penaltyAreaInit;
      for(uint PenaltyIter = 0; PenaltyIter < 100; PenaltyIter++)
      {
	double RegEnergy = 1.0*kvisc;
	uint   RegIter = 0;
	double RegTol = viscTol*(shellBody->totalStrainEnergy());
	
	while (RegEnergy > RegTol && RegIter < 100)
	{
	  solver.setBounds(BoundType, LowerBound, UpperBound);
	  solver.solve(&model);
	  RegEnergy = RegularizationBody.totalStrainEnergy();
	  cout << "              PenaltyIter-RegIter = " << 
	    PenaltyIter << " - " << RegIter  << " RegEnergy = " << RegEnergy << endl;
	  RegularizationBody.resetRefConf();
	  RegIter++;
	}

	// Update Lagrange multipliers
	shellBody->updateFixedPressure();
	shellBody->updateFixedTension();
	// Update Penalty to improve convergence
	shellBody->updatePenaltyVolumeArea(penaltyVolume*=3, penaltyArea*=3, 0.0);
	
	double vred = 6.0*sqrt(M_PI)*shellBody->volume()/pow(shellBody->area(),3.0/2.0);

	cout << "       PenaltyIter               = " << PenaltyIter << ":"    << endl
	     << "       Computed   reduced volume = " << vred                  << endl
	     << "       Prescribed reduced volume = " << Vred                  << endl
	     << "       Current pressure          = " << shellBody->pressure() << endl
	     << "       Fixed pressure            = " << shellBody->fixedPressure() << endl
	     << "       Current tension           = " << shellBody->tension()  << endl
	     << "       Fixed tension             = " << shellBody->fixedTension()  << endl
	     << "       RegTol                    = " << RegTol                 << endl;
		
	if( (abs(vred-Vred) < penaltyVolumeTol || Vmin > 1.0) && abs(FixedShellArea - shellBody->area()) < penaltyAreaTol )
	{
	  break;
	}
	
      }
      
      Dislocation = PrintVirus.printMaster(++PrintCount);
      sprintf(vtkname,"Capsid-%04d",PrintCount);
      shellBody->printParaview(vtkname);

      // Compute Force vs Displacement
      Vector3D Force0(0.0), Force1(0.0);
      for (uint i = 0; i < defNodesV0.size(); i++) {
	Force0 += defNodesV0[i]->force();
      }
      for (uint i = 0; i < defNodesV1.size(); i++) {
	Force1 += defNodesV1[i]->force();
      }
      cout << "     u = " << defNodesV0[0]->getPoint(2) << "   "
	   << "Force0 = " << Force0(0) << " " << Force0(1) << " " <<Force0(2) << "   "
	   << "Force1 = " << Force1(0) << " " << Force1(1) << " " <<Force1(2) << endl;
      cout << "Pressure = " << shellBody->pressure() << " Surface Tension = " << shellBody->tension() << endl;
      cout << "Vertices valence [5, 6, 7] = " << Dislocation[0] << " " << Dislocation[1] << " " << Dislocation[2] << endl;

    }
  } // end of restart (if any)
  

  cout << " -------------- Volume has been reduced -------------- " << endl;
  cout << " Shell body volume = " << shellBody->volume() << endl;
  cout << " Shell body area   = " << shellBody->area()   << endl;
  cout << " --------------------- ----------------------- ------- " << endl;

  


  
  
  
  


  // Apply force by controlled displacements
  double DispStep = MaxDisp/nDispSteps;
  Vred = Vmin;
  shellBody->reduceVolume(Vred); // It should not be necessary if restart is not used since _prescribed volume is not changed
  cout << endl << " Average edge length = " << shellBody->AverageEdgeLength() << endl << endl;

  // Erase BC used during volume reduction
  BoundType = 0;
  LowerBound = 0.0;
  UpperBound = 0.0;

  for (uint step = 0; step < nDispSteps; step++) {

    DeformationNode<3>::Point x;
    NodeBase::DofIndexMap idx(3);
    // Apply displacement
    for (uint i = 0; i < defNodesV0.size(); i++) {
      x = defNodesV0[i]->point();
      idx = defNodesV0[i]->index();
      // Constrained applied in X, Y, and Z directions
	BoundType(idx[0]) = 2;
	LowerBound(idx[0]) = x(0);
	UpperBound(idx[0]) = x(0);

	BoundType(idx[1]) = 2;
	LowerBound(idx[1]) = x(1);
	UpperBound(idx[1]) = x(1);

	BoundType(idx[2]) = 2;
	LowerBound(idx[2]) = x(2) + DispStep;
	UpperBound(idx[2]) = x(2) + DispStep;
    }
    for (uint i = 0; i < defNodesV1.size(); i++) {
      x = defNodesV1[i]->point();
      idx = defNodesV1[i]->index();
      // Constrained applied in X, Y, and Z directions
        BoundType(idx[0]) = 2;
        LowerBound(idx[0]) = x(0);
        UpperBound(idx[0]) = x(0);
      
        BoundType(idx[1]) = 2;
        LowerBound(idx[1]) = x(1);
        UpperBound(idx[1]) = x(1);
      
	BoundType(idx[2]) = 2;
	LowerBound(idx[2]) = x(2) - DispStep;
	UpperBound(idx[2]) = x(2) - DispStep;
    }

    double penaltyVolume = penaltyVolumeInit;
    double penaltyArea   = penaltyAreaInit;

    // Update neighbors used in ProteinBody
    if (epsilon > TolToAddProteinBody) {
      ProteinBody->recomputeNeighbors(ScaledSearchR);
    }
    cout << endl << "              Average edge length = " << shellBody->AverageEdgeLength() << endl << endl;
    for(uint PenaltyIter = 0; PenaltyIter < 100; PenaltyIter++)
    {
      double RegEnergy = 1.0*kvisc;
      uint   RegIter = 0;
      double RegTol = viscTol*(shellBody->totalStrainEnergy() );
      // + fabs(ProteinBody.totalStrainEnergy()) );
      while (RegEnergy > RegTol && RegIter < 100)
      {
	solver.setBounds(BoundType, LowerBound, UpperBound);
	solver.solve(&model);
	RegEnergy = RegularizationBody.totalStrainEnergy();
	cout << "              PenaltyIter-RegIter = " << 
	  PenaltyIter << " - " << RegIter  << " RegEnergy = " << RegEnergy << 
	  " ShellBodyEnergy = " << shellBody->totalStrainEnergy() <<
          " ProteinBodyEnergy = " << ProteinBody->totalStrainEnergy() << endl;
	RegularizationBody.resetRefConf();
	RegIter++;
      } // Viscous regularization loop

      // Update Lagrange multipliers
      shellBody->updateFixedPressure();
      shellBody->updateFixedTension();
      // Update Penalty to improve convergence
      shellBody->updatePenaltyVolumeArea(penaltyVolume*=3, penaltyArea*=3, 0.0);
      
      double vred = 6.0*sqrt(M_PI)*shellBody->volume()/pow(shellBody->area(),3.0/2.0);

      cout << "       PenaltyIter               = " << PenaltyIter << ":"    << endl
	   << "       Computed   reduced volume = " << vred                  << endl
	   << "       Prescribed reduced volume = " << Vred                  << endl
	   << "       Current pressure          = " << shellBody->pressure() << endl
	   << "       Fixed pressure            = " << shellBody->fixedPressure() << endl
	   << "       Current tension           = " << shellBody->tension()  << endl
	   << "       Fixed tension             = " << shellBody->fixedTension()  << endl;
		
      if( (abs(vred-Vred) < penaltyVolumeTol || Vmin > 1.0) && abs(FixedShellArea - shellBody->area()) < penaltyAreaTol )
      {
	break;
      }

    } // Penalty enforcement loop

    // Print displacement loop results
    Dislocation = PrintVirus.printMaster(++PrintCount);
    sprintf(vtkname,"Capsid-%04d",PrintCount);
    shellBody->printParaview(vtkname);

    // Compute Force vs Displacement
    Vector3D Force0(0.0), Force1(0.0);
    for (uint i = 0; i < defNodesV0.size(); i++) {
      Force0 += defNodesV0[i]->force();
    }
    for (uint i = 0; i < defNodesV1.size(); i++) {
      Force1 += defNodesV1[i]->force();
    }
    cout << "     u = " << DispStep*(step+1) << "   "
         << "Force0 = " << Force0(0) << " " << Force0(1) << " " <<Force0(2) << "   "
	 << "Force1 = " << Force1(0) << " " << Force1(1) << " " <<Force1(2) << endl;
    cout << "Pressure = " << shellBody->pressure() << " Surface Tension = " << shellBody->tension() << endl;
    cout << "Vertices valence [5, 6, 7] = " << Dislocation[0] << " " << Dislocation[1] << " " << Dislocation[2] << endl;



    // Remesh if necessary
    uint ElementsChanged = shellBody->Remesh(ARtol, bodyMaterial, quadOrder);
    if (ElementsChanged > 0) {
      // // Print remeshed - Force is measured before remeshing
      // PrintVirus.printMaster(++PrintCount);
      // sprintf(vtkname,"Capsid-%04d",PrintCount);
      // shellBody->printParaview(vtkname);
      cout << "     ElementsChanged = " << ElementsChanged << endl;
    }



  } // Controlled displacement loop


   

  
     
  
  





  // -------------------------------------------------------------------
  // Output
  // -------------------------------------------------------------------
  // Calculate capsid center, average radius and asphericity
  tvmet::Vector<double,3> center(0.0);
  Ravg = 0.0;
  double deltaR = 0.0, Rtemp = 0.0, asphericity = 0.0;
  for(uint i = 0; i < defNodes.size(); i++)
  {
    center += defNodes[i]->point();
    Ravg   += tvmet::norm2(defNodes[i]->point());
  }
  center /= defNodes.size();
  Ravg /= defNodes.size();
  cout << endl << "Center      = " << center  << endl;
  
  for(uint i = 0; i < defNodes.size(); i++)
  {
    Rtemp = tvmet::norm2(defNodes[i]->point());
    deltaR += pow(Rtemp-Ravg,2);
  }
  asphericity = sqrt( deltaR/defNodes.size() )/Ravg;
   

  
  //---------------//
  // Print results //
  // Print the total energies of three bodies
  cout << " WshellBody                 = " << shellBody->totalStrainEnergy() << endl;
  // cout << " WpotentialBody_1Ring= " << potentialBody_FirstRing.totalStrainEnergy() << endl;
  // cout << " WpotentialBody_2Ring= " << potentialBody_SecondRing.totalStrainEnergy() << endl;
  // cout << " Total Energy        = " << (shellBody->totalStrainEnergy() + potentialBody_FirstRing.totalStrainEnergy()) << endl;
 //+ potentialBody_SecondRing.totalStrainEnergy()) << endl;
  cout << " Regularization body energy = " << RegularizationBody.totalStrainEnergy() << endl;
  cout << " Final shell body volume    = " << shellBody->volume() << endl;
  cout << " Final shell body area      = " << shellBody->area()   << endl;
  cout << "---------------------------------------------------"   << endl;
  
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
  





  
  // -------------------------------------------------------------------
  // Clean up
  // -------------------------------------------------------------------
delete shellBody;
delete ProteinBody;
  
  for (uint i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }
  
  return (0);  
}
