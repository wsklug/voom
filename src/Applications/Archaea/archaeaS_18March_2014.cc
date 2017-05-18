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
#include "LennardBody.h"
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
  int restart = false;
  double SphereR = 1.0;
 
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

  // Penalty parameters
  double penaltyVolumeInit = 1.0e3;
  double penaltyAreaInit   = 1.0e3;

  // Potential input parameters
  double PotentialSearchRF = 1.0;
  double epsilon = 1.0;
  double sigma = 1.0;

  // Connectivity search
  uint RconnSF = 1.0;

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
  inp >> temp >> SphereR;
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> ForceSearchR;
  inp >> temp >> MaxDisp;
  inp >> temp >> nDispSteps;
  inp >> temp >> kvisc;
  inp >> temp >> penaltyVolumeInit;
  inp >> temp >> penaltyAreaInit;
  inp >> temp >> PotentialSearchRF;
  inp >> temp >> epsilon;
  inp >> temp >> sigma;
  inp >> temp >> RconnSF;

  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " outputFileName          : " << outputFileName    << endl
       << " restart                 : " << restart           << endl
       << " SphereR                 : " << SphereR           << endl
       << " KC                      : " << KC                << endl
       << " KG                      : " << KG                << endl
       << " C0                      : " << C0                << endl
       << " Force Search Radius     : " << ForceSearchR      << endl
       << " Max Disp                : " << MaxDisp           << endl
       << " N disp steps            : " << nDispSteps        << endl
       << " Viscous regularization  : " << kvisc             << endl
       << " Penalty volume factor   : " << penaltyVolumeInit << endl
       << " Penalty area factor     : " << penaltyAreaInit   << endl  
       << " Potential Search factor : " << PotentialSearchRF << endl
       << " epsilon                 : " << epsilon           << endl
       << " sigma                   : " << sigma             << endl
       << " Rconn scaling factr     : " << RconnSF           << endl;
      
  






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
  vector<DeformationNode<3>* > defNodes, defNodesV0, defNodesV1;
  double Ravg = 0;

  // Opposite vertices where forces are applied
  Vector3D V0(0.0, 0.0, 1.0);
  Vector3D V1(0.0, 0.0,-1.0);
  uint CountV0nodes = 0, CountV1nodes = 0, CountSolverNodes = 0;

  // Bounds for solver
    blitz::Array<double,1> LowerBound, UpperBound;
    blitz::Array<int,1> BoundType;

   
  
  // Input .vtk or inp files containing nodes and connectivities
  string token;
  if (restart) {
    ifs >> token;
    while( token != "POINTS" ) ifs >> token;
	cout << token << endl;
  }
  else {
    ifs >> token;
  }

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

  if (restart) {
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
    }
    Ravg /= nodes.size();
  }
  else {
    // read in points
    for(uint i = 0; i < npts; i++) {
      uint NodeNum = 0;
      DeformationNode<3>::Point x;
      ifs >> NodeNum >> x(0) >> x(1) >> x(2);
      Ravg += tvmet::norm2(x);
      
      NodeBase::DofIndexMap idx(3);
      for(uint j = 0; j < 3; j++) { 
	idx[j] = dof++;
      }
      DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
      
      nodes.push_back( n );
      defNodes.push_back( n );
      
      // Distinguish between free and prescribed nodes
      if (tvmet::norm2(x-V0) < ForceSearchR) {
	defNodesV0.push_back(n);
      }
      else if (tvmet::norm2(x-V1) < ForceSearchR) {
	defNodesV1.push_back(n);
      }
    } // end of reading points
  }
  CountV0nodes = defNodesV0.size();
  CountV1nodes = defNodesV1.size();
  
  Ravg /= nodes.size();
  cout << endl << "V0 nodes number = " << CountV0nodes << endl;
  cout << endl << "V1 nodes number = " << CountV1nodes << endl;
  // assert(CountV0nodes == CountV1nodes);
  

   
  // Read in triangle connectivities
  
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  uint ntri = 0, ElemNum = 0, tmp = 0;
  double AverageEdgeLength = 0.0;

  if(restart) {
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
  }
  else {
    ifs >> token;
    ifs >> ntri;
    connectivities.reserve(ntri);
    for (uint i = 0; i < ntri; i++)
    {
      ifs >> ElemNum;
      ifs >> ct(0) >> ct(1) >> ct(2);
      ct(0) -= 1;  ct(1) -= 1;  ct(2) -= 1;
      // Compute initial average edge length - used in potential body
      AverageEdgeLength += tvmet::norm2(defNodes[ct(0)]->point()- defNodes[ct(1)]->point());
      AverageEdgeLength += tvmet::norm2(defNodes[ct(1)]->point()- defNodes[ct(2)]->point());
      AverageEdgeLength += tvmet::norm2(defNodes[ct(2)]->point()- defNodes[ct(0)]->point());

      connectivities.push_back(ct);
    }
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
  int quadOrder = 2;
  
  // Bending material object, to be copied when body generates new elements
  EvansElastic bodyMaterial(KC, KG, C0, 0.0, 0.0);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  double pressure = 0.0;//-10.0;
  double tension  = 0.0;//-5.0;
  /* LoopShellBody<EvansElastic> * shellBody = new LoopShellBody<EvansElastic>(bodyMaterial, 
									    connectivities, 
									    nodes, 
									    quadOrder);
  */
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
  

  // Initialize body
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  cout << "Initial shell body volume = " << shellBody->volume() << endl;
  cout << "Initial shell body area   = " << shellBody->area()   << endl;
  // shellBody->checkConsistency();
  // Set area body to the one of a perfect sphere with radius R
  double PerfectSphereArea = 4.0*M_PI*pow(SphereR, 2.0);
  shellBody->setPrescribedArea(PerfectSphereArea);


  // Create viscosity body
  ViscosityBody RegularizationBody(defNodes, connectivities, kvisc);
  RegularizationBody.compute(true, false, false);
  cout << "Initial regularization body energy = " << RegularizationBody.totalStrainEnergy() << endl;
  cout << "Initial regularization body energy = " << RegularizationBody.energy() << endl;




  

  // Initiliaze Lennard-Jones potential bodies
  
  double ScaledSearchR = PotentialSearchRF*AverageEdgeLength;
  double ScaledSigma   = sigma*AverageEdgeLength;
  LennardBody LennardPotentialBody(defNodes, ScaledSearchR, epsilon, ScaledSigma);
  // potentialBody->checkConsistency();
  
  LennardPotentialBody.compute(true, false, false);
  cout << "Initial potential body (FirstRing) energy = " << LennardPotentialBody.totalStrainEnergy() << endl;
  cout << "Initial potential body (FirstRing) energy = " << LennardPotentialBody.energy() << endl;
  
  
  /*
  // First ring interactions
  double ScaledSearchR = PotentialSearchRF*AverageEdgeLength;
  double ScaledSigma   = sigma*AverageEdgeLength;
  LennardBody potentialBody_FirstRing(defNodes, ScaledSearchR, epsilon, ScaledSigma);
  // potentialBody->checkConsistency();
  
  potentialBody_FirstRing.compute(true, false, false);
  cout << "Initial potential body (FirstRing) energy = " << potentialBody_FirstRing.totalStrainEnergy() << endl;
  cout << "Initial potential body (FirstRing) energy = " << potentialBody_FirstRing.energy() << endl;

  // Second ring interactions
  LennardBody potentialBody_SecondRing(defNodes, 2.0*ScaledSearchR, 0.1*epsilon, ScaledSigma);
  // potentialBody->checkConsistency();
  
  potentialBody_SecondRing.compute(true, false, false);
  cout << "Initial potential body (SecondRing) energy = " << potentialBody_SecondRing.totalStrainEnergy() << endl;
  cout << "Initial potential body (SecondRing) energy = " << potentialBody_SecondRing.energy() << endl;
  */



  // Create Model and solver
  Model::BodyContainer bdc;
  bdc.push_back(shellBody);
  // bdc.push_back(LennardPotentialBody);
  bdc.push_back(&RegularizationBody);
  // bdc.push_back(&potentialBody_FirstRing);
  // bdc.push_back(&potentialBody_SecondRing);
  
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
  PrintingArchaeaS PrintVirus(modelName, outputFileName+"iter", bdc, double(RconnSF*AverageEdgeLength), defNodes);

  // Print before everything
  uint PrintCount = -1;
  PrintVirus.printMaster(PrintCount);
  char vtkname[100]; 
  sprintf(vtkname,"Capsid-%04d",PrintCount);
  shellBody->printParaview(vtkname);







  DeformationNode<3>::Point x;
  for (uint i = 0; i<defNodes.size(); i++) {
  // Introduce bias in initial configuration
    x = defNodes[i]->point();
    x(0) *= (0.9 + fabs(x(2))*0.1);
    x(1) *= (0.9 + fabs(x(2))*0.1);
    x(2) *= (1.0 + fabs(x(2))*0.1);
    defNodes[i]->setPoint(x);
  }
  RegularizationBody.resetRefConf();






  // -------------------------------------------------------------------
  // Solve
  // -------------------------------------------------------------------

  // Set bounds on nodes from which to pull
  NodeBase::DofIndexMap idx(3);
  // Apply constraint
  for (uint i = 0; i < defNodesV0.size(); i++) {
    x = defNodesV0[i]->point();
    idx = defNodesV0[i]->index();
    // Constrained applied in x, y, z directions
    BoundType(idx[0]) = 2;
    LowerBound(idx[0]) = x(0);
    UpperBound(idx[0]) = x(0);

    BoundType(idx[1]) = 2;
    LowerBound(idx[1]) = x(1);
    UpperBound(idx[1]) = x(1);

    BoundType(idx[2]) = 2;
    LowerBound(idx[2]) = x(2);
    UpperBound(idx[2]) = x(2);
  }
  for (uint i = 0; i < defNodesV1.size(); i++) {
    x = defNodesV1[i]->point();
    idx = defNodesV1[i]->index();

    // Initially the constrained is applied only in the X and Y directions
    BoundType(idx[0]) = 2;
    LowerBound(idx[0]) = x(0);
    UpperBound(idx[0]) = x(0);

    BoundType(idx[1]) = 2;
    LowerBound(idx[1]) = x(1);
    UpperBound(idx[1]) = x(1);
  }

  // Solve to minimize the energy - Initial relaxed configuration
  // solver.setBounds(BoundType, LowerBound, UpperBound);
  // solver.solve(&model);

  PrintVirus.printMaster(++PrintCount);
  sprintf(vtkname,"Capsid-%04d",PrintCount);
  shellBody->printParaview(vtkname);



  // Solve to set fixed volume and fixed area
  double PenaltyTol = 1.0e-3;
  // Initiliaze shell by reducing its volume to 0.94
  double Vred = 1.0;

  if (!restart)
  {
    for(Vred = 0.99; Vred >= 0.9; Vred -= 0.01)
    {
      shellBody->reduceVolume(Vred);
      
      double penaltyVolume = penaltyVolumeInit;
      double penaltyArea   = penaltyAreaInit;
      for(uint PenaltyIter = 0; PenaltyIter < 10; PenaltyIter++)
      {
	double RegEnergy = 1.0*kvisc;
	uint   RegIter = 0;
	double RegTol = 1.0e-3*shellBody->totalStrainEnergy();
	while (RegEnergy > RegTol && RegIter < 10)
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
	     << "       Fixed tension             = " << shellBody->fixedTension()  << endl;
		
	if( abs(vred-Vred) < PenaltyTol )
	{
	  break;
	}
	
      }
      PrintVirus.printMaster(++PrintCount);

      sprintf(vtkname,"Capsid-%04d",PrintCount);
      shellBody->printParaview(vtkname);
    }
  } // end of restart (if any)
  

  cout << " -------------- Volume has been reduced -------------- " << endl;
  cout << " Shell body volume = " << shellBody->volume() << endl;
  cout << " Shell body area   = " << shellBody->area()   << endl;
  cout << " --------------------- ----------------------- ------- " << endl;

  


  
  // Apply force by controlled displacements
  double DispStep = MaxDisp/nDispSteps;
  Vred = 0.9;

  for (uint step = 0; step < nDispSteps; step++) {
    // shellBody->reduceVolume(Vred); // Not necessary since _prescribed volume is not changed

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
    for(uint PenaltyIter = 0; PenaltyIter < 10; PenaltyIter++)
    {
      double RegEnergy = 1.0*kvisc;
      uint   RegIter = 0;
      double RegTol = 1.0e-3*shellBody->totalStrainEnergy();
      while (RegEnergy > RegTol && RegIter < 10)
      {
	solver.setBounds(BoundType, LowerBound, UpperBound);
	solver.solve(&model);
	RegEnergy = RegularizationBody.totalStrainEnergy();
	cout << "              PenaltyIter-RegIter = " << 
	  PenaltyIter << " - " << RegIter  << " RegEnergy = " << RegEnergy << endl;
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
		
      if( abs(vred-Vred) < PenaltyTol )
      {
	break;
      }

    } // Penalty enforcement loop

    // Print displacement loop results
    PrintVirus.printMaster(++PrintCount);
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
  
  for (uint i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }
  
  return (0);  
}
