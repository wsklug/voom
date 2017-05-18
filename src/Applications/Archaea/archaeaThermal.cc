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

#include "BrownianDynamics3D.h"

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
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;

  // Applied force parameters
  double ForceSearchR = 0.0;

  // Penalty parameters
  double penaltyVolumeInit = 1.0e3;
  double penaltyAreaInit   = 1.0e3;
  double penaltyVolumeTol  = 1.0e-3;
  double penaltyAreaTol    = 1.0e-3;
  double InitialPressure   = 0.0;
  double InitialTension    = 0.0;

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
  double Vmin = 0.93;

  // PotentialBodyFull flag
  uint PotentialBodyFull = 1;  // 0 -> nodes at which we apply displacements are not proteins (no Morse potential interactions)
                               // 1 -> nodes at which we apply displacements are proteins (Included in potential body)

  int quadOrder = 1;

  double mobility = 1.0;

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
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> ForceSearchR;
  inp >> temp >> penaltyVolumeInit;
  inp >> temp >> penaltyAreaInit;
  inp >> temp >> penaltyVolumeTol;
  inp >> temp >> penaltyAreaTol;
  inp >> temp >> InitialPressure;
  inp >> temp >> InitialTension;
  inp >> temp >> PotentialSearchRF;
  inp >> temp >> epsilon;
  inp >> temp >> sigma;
  inp >> temp >> Rshift;
  inp >> temp >> RconnSF;
  inp >> temp >> ARtol;
  inp >> temp >> Vmin;
  inp >> temp >> PotentialBodyFull;
  inp >> temp >> quadOrder;
  inp >> temp >> mobility;


  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " outputFileName          : " << outputFileName    << endl
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
       << " Initial pressure        : " << InitialPressure   << endl
       << " Initial tension         : " << InitialTension    << endl  
       << " Potential Search factor : " << PotentialSearchRF << endl
       << " epsilon                 : " << epsilon           << endl
       << " sigma                   : " << sigma             << endl
       << " Rshift                  : " << Rshift            << endl
       << " Rconn scaling factr     : " << RconnSF           << endl
       << " Aspect Ratio tol        : " << ARtol             << endl
       << " Min reduced volume      : " << Vmin              << endl
       << " PotentialBodyFull       : " << PotentialBodyFull << endl
       << " Loop Shell quadOrder    : " << quadOrder         << endl
       << " mobility                : " << mobility          << endl;


      
  






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
  vector<BrownianNode<3>* > defNodes, defNodesV0, defNodesV1, SolverNodes;
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


  ifs >> token;   // skip number type in vtk file
  // read in points
  for(uint i = 0; i < npts; i++) {
    BrownianNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) idx[j] = dof++;
    
    BrownianNode<3>* n = new BrownianNode<3>(i,idx,x);
    BrownianNode<3>::Matrix MobilityMatrix(0.0);
    MobilityMatrix(0,0) = mobility;
    MobilityMatrix(1,1) = mobility;
    MobilityMatrix(2,2) = mobility;
    n->setMobility(MobilityMatrix);
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
    else {
      SolverNodes.push_back( n );
    }
   
  }
  Ravg /= nodes.size();
  
  CountV0nodes = defNodesV0.size();
  CountV1nodes = defNodesV1.size();
  
  cout << endl << "V0 nodes number = " << CountV0nodes << endl;
  cout << endl << "V1 nodes number = " << CountV1nodes << endl;
  

   
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
    
    connectivities.push_back(ct);
  }
 
  
  
  cout << "Number of nodes     = " << nodes.size()          << endl
       << "Number of elements  = " << connectivities.size() << endl
       << "Ravg                = " << Ravg                  << endl;

  // Close mesh file in vtk or inp format
  ifs.close();
  


  // -------------------------------------------------------------------
  // Create model
  // -------------------------------------------------------------------
  // Create Loop shell body with both bending and in-plane elastic energy
  
  // Bending material object, to be copied when body generates new elements
  EvansElastic bodyMaterial(KC, KG, C0, 0.0, 0.0);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  double pressure = InitialPressure;
  double tension  = InitialTension;
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
									    augmented, 
									    augmented,
									    noConstraint);
  shellBody->reduceVolume(Vmin); 

  // Initialize body
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  cout << "Initial shell body volume = " << shellBody->volume() << endl;
  cout << "Initial shell body area   = " << shellBody->area()   << endl;
  cout << "Average edge length from body = " << shellBody->AverageEdgeLength() << endl << endl;
  // shellBody->checkConsistency();
  


  // Initiliaze potential body
	Morse ProteinPotential(epsilon, sigma, Rshift); 

	// Then initialize potential body
	PotentialBody * ProteinBody;
	if (PotentialBodyFull == 1) {
	  ProteinBody = new PotentialBody(&ProteinPotential, defNodes, PotentialSearchRF); }
	else {
	  ProteinBody = new PotentialBody(&ProteinPotential, SolverNodes, PotentialSearchRF);
	}
	  
  	// potentialBody->checkConsistency();
  
  	ProteinBody->compute(true, false, false);
  	cout << "Initial protein body (FirstRing) energy = " << ProteinBody->totalStrainEnergy() << endl;
  	cout << "Initial protein body (FirstRing) energy = " << ProteinBody->energy() << endl;







  
  // Set solver and solve
  int printStride =  -1;
  int nSteps = 100;
  int nPrintSteps = 10;
  double dt = 1.0e-1;

  BrownianDynamics3D * brownian = new BrownianDynamics3D(SolverNodes, printStride);

  brownian->pushBackBody(shellBody);
  Model::BodyContainer bdc;
  bdc.push_back(shellBody);
  double TolToAddProteinBody = 1.0e-12;
  if (epsilon > TolToAddProteinBody) {
    brownian->pushBackBody(ProteinBody);
    bdc.push_back(ProteinBody);
  }
  else {
    cout << "No protein body added." << endl;
  }







  // Initialize printing utils
  cout << "Seraching connectivity over R = " << RconnSF << endl;
  PrintingArchaeaS PrintVirus(modelName, outputFileName+"iter", bdc, RconnSF, defNodes, &ProteinPotential);

  // Print before everything
  uint PrintCount = 0;
  vector<uint > Dislocation(3,0);
  Dislocation = PrintVirus.printMaster(PrintCount);
  char vtkname[100]; 
  sprintf(vtkname,"Capsid-%04d",PrintCount);
  shellBody->printParaview(vtkname);

  






  // -------------------------------------------------------------------
  // Solve
  // -------------------------------------------------------------------

  for(int printStep=0; printStep<=nPrintSteps; printStep++) {
     
	double penaltyVolume = penaltyVolumeInit;
    	double penaltyArea   = penaltyAreaInit;

	brownian->run( nSteps, dt );

	// Update neighbors used in ProteinBody
    	if (epsilon > TolToAddProteinBody) {
      		ProteinBody->recomputeNeighbors(PotentialSearchRF);
    	}
    	cout << endl << "              Average edge length = " << shellBody->AverageEdgeLength() << endl << endl;


    	for(uint PenaltyIter = 0; PenaltyIter < 100; PenaltyIter++)
    	{
		// Update Lagrange multipliers
      		shellBody->updateFixedPressure();
      		shellBody->updateFixedTension();
      		// Update Penalty to improve convergence
      		shellBody->updatePenaltyVolumeArea(penaltyVolume*=3, penaltyArea*=3, 0.0);
      
      		double vred = 6.0*sqrt(M_PI)*shellBody->volume()/pow(shellBody->area(),3.0/2.0);

		brownian->run( nSteps, dt );

      		cout << "       PenaltyIter               = " << PenaltyIter << ":"    << endl
	   	<< "       Computed   reduced volume = " << vred                  << endl
	   	<< "       Prescribed reduced volume = " << Vmin                  << endl
	   	<< "       Current pressure          = " << shellBody->pressure() << endl
	   	<< "       Fixed pressure            = " << shellBody->fixedPressure() << endl
	   	<< "       Current tension           = " << shellBody->tension()  << endl
	   	<< "       Fixed tension             = " << shellBody->fixedTension()  << endl;
		
      		if( (abs(vred-Vmin) < penaltyVolumeTol || Vmin > 1.0) && abs(FixedShellArea - shellBody->area()) < penaltyAreaTol )
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
    	cout << "     u = " << printStep << "   "
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
}






    std::cout << "BrownianDynamics: step "<< printStep*nSteps
              << std::setprecision( 16 )
              << " | energy = "
              << brownian->energy() << std::endl;



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
