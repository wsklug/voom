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
#include "LoopShell.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "VoomMath.h"
#include "Utils/PrintingArchaea.h"

using namespace voom;
using namespace std;

int main(int argc, char* argv[])
{
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
  string OutlineName;
  
  // Minimization constants
  int min = 0;            // 0 -> no force, no Leonard-Johns potential, no edge swapping (1 -> force, 2 -> force+Leonard-Johns, 3 -> force+Leonard=Johns+edge-swapping)
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  // Stretching
  double mu = 1.0;
  double kS = 1.0;
  // Applied force parameters
  double ForceSearchR = 0.0;
  double MaxDisp = 1.0;
  double nDispSteps = 1.0;
  // Potential input parameters
  double PotentialSearchRF = 1.0;
  double epsilon = 1.0;
  double sigma = 1.0;


  // Output
  unsigned int NPout = 0; // Nodes needed to build the outline
  string OutFile;
  int refinement = 1;

  // Reading input from file passed as argument
  ifstream inp;
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  string temp;
  
  inp >> temp >> modelName;
  inp >> temp >> OutlineName;
  inp >> temp >> min;
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> mu;
  inp >> temp >> kS;
  inp >> temp >> NPout;
  inp >> temp >> OutFile;
  inp >> temp >> refinement;
  inp >> temp >> ForceSearchR;
  inp >> temp >> MaxDisp;
  inp >> temp >> nDispSteps;
  inp >> temp >> PotentialSearchRF;
  inp >> temp >> epsilon;
  inp >> temp >> sigma;

  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " OutlineName             : " << OutlineName       << endl
       << " min                     : " << min               << endl
       << " KC                      : " << KC                << endl
       << " KG                      : " << KG                << endl
       << " C0                      : " << C0                << endl
       << " mu                      : " << mu                << endl
       << " kS                      : " << kS                << endl
       << " NPout                   : " << NPout             << endl
       << " OutFile                 : " << OutFile           << endl
       << " Refinement              : " << refinement        << endl
       << " Force Search Radius     : " << ForceSearchR      << endl
       << " Max Disp                : " << MaxDisp           << endl
       << " N disp steps            : " << nDispSteps        << endl
       << " Potential Search factor : " << PotentialSearchRF << endl
       << " epsilon                 : " << epsilon           << endl
       << " sigma                   : " << sigma             << endl;
      
  





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
  double phi = 0.5*(1 + sqrt(5.0));
  Vector3D V0(0.0, 1.0,-phi);
  Vector3D V1(0.0,-1.0, phi);
  uint CountV0nodes = 0, CountV1nodes = 0, CountSolverNodes = 0;

  // Bounds for solver
    blitz::Array<double,1> LowerBound, UpperBound;
    blitz::Array<int,1> BoundType;

   
  
  // Input .vtk file containing nodes and connectivities
  string token;
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
     BoundType  = 0;  // Default is no bound
  ifs >> token;   // skip number type

  // read in points
  for(uint i = 0; i < npts; i++) {
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
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
      /*for(uint j = 0; j < 3; j++) {
	BoundType(idx[j]) = 2; // Both upper and lower bound
	LowerBound(idx[j]) = x(j);
	UpperBound(idx[j]) = x(j);
	}*/
      defNodesV0.push_back(n);
    }
    else if (tvmet::norm2(x-V1) < ForceSearchR) {
      /*for(uint j = 0; j < 3; j++) {
	BoundType(idx[j]) = 2; // Both upper and lower bound
	LowerBound(idx[j]) = x(j);
	UpperBound(idx[j]) = x(j);
	}*/
      defNodesV1.push_back(n);
    }
  } // end of reading points
  CountV0nodes = defNodesV0.size();
  CountV1nodes = defNodesV1.size();
  
  Ravg /= nodes.size();
  cout << endl << "V0 nodes number = " << CountV0nodes << endl;
  cout << endl << "V1 nodes number = " << CountV1nodes << endl;
  assert(CountV0nodes == CountV1nodes);
  

   
  // Read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  int ntri = 0, tmp = 0;
  double AverageEdgeLength = 0.0;
  ifs >> ntri;
  connectivities.reserve(ntri);
  
  ifs >> tmp;
  for (uint i = 0; i < ntri; i++)
  {
    ifs >> tmp;
    if(tmp != 3) { 
      cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    }
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
  
  // Close mesh file in vtk format
  ifs.close();
  
 



  

  // Create Loop shell body with both bending and in-plane elastic energy
  int quadOrder = 1;
  
  // Bending material object, to be copied when body generates new elements
  EvansElastic bodyMaterial(KC, KG, C0, mu, kS);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  LoopShellBody<EvansElastic> * shellBody = new LoopShellBody<EvansElastic>(bodyMaterial, connectivities, nodes, quadOrder);
  // Initialize body
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  // shellBody->checkConsistency();

  // Initiliaze Lennard-Jones potential bodies
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




  // Create model and set solver
  switch (min) {
  case 1: // no force - no Leonnard-Johns potential - no edge-swapping
    // Nothing to be done, all nodes are DOF
    break;
  case 2: // force - no Leonnard-Johns potential - no edge-swapping
    // TODO
    break;
  }

  // Create Model
  Model::BodyContainer bdc;
  bdc.push_back(shellBody);
  bdc.push_back(&potentialBody_FirstRing);
  bdc.push_back(&potentialBody_SecondRing);
  
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
  // Set bounds on nodes from which to pull
  // solver.setBounds(BoundType, LowerBound, UpperBound);

  cout << "Solver DOF = " << model.dof() << endl;

  // Compute FvK number before solving
  double FvK = 0.0, Y = 0.0;
  Y = 4.0*kS*mu/(kS+mu);
  FvK = Y*Ravg*Ravg/KC;
  cout << " FvK number before solving = " << FvK           << endl;

  // Solve to minimize the energy - Initial relaxed configuration
  solver.solve(&model);
  PrintingArchaea PrintVirus(modelName, OutlineName, NPout, OutFile+"iter", OutFile+"outline", bdc, refinement);
  PrintVirus.printMaster(0);

  // Apply force by controlled displacements
  double DispStep = MaxDisp/nDispSteps;
  Vector3D V01(0.0);
  V01 = V1-V0; V01 = V01/tvmet::norm2(V01);
  for (uint step = 0; step < nDispSteps; step++) {
    DeformationNode<3>::Point x;
    NodeBase::DofIndexMap idx(3);
    // Apply displacement
    for (uint i = 0; i < defNodesV0.size(); i++) {
      x = defNodesV0[i]->point();
      idx = defNodesV0[i]->index();
      for(uint j = 0; j < 3; j++) {
	BoundType(idx[j]) = 2;
	LowerBound(idx[j]) = x(j) - V01(j)*DispStep;
	UpperBound(idx[j]) = x(j) - V01(j)*DispStep;;
      }
    }
    for (uint i = 0; i < defNodesV1.size(); i++) {
      x = defNodesV1[i]->point();
      idx = defNodesV1[i]->index();
      for(uint j = 0; j < 3; j++) {
	BoundType(idx[j]) = 2;
	LowerBound(idx[j]) = x(j) + V01(j)*DispStep;
	UpperBound(idx[j]) = x(j) + V01(j)*DispStep;;
      }
    }
    solver.setBounds(BoundType, LowerBound, UpperBound);
    solver.solve(&model);
    PrintVirus.printMaster(step);
  }
  
  





  // Output
  // Calculate capsid center, average radius and asphericity
  tvmet::Vector<double,3> center(0.0);
  Ravg = 0.0;
  double deltaR = 0.0, Rtemp = 0.0, asphericity = 0.0;
  for(uint i = 0; i < defNodes.size(); i++)
  {
    center += defNodes[i]->point();
    Ravg += tvmet::norm2(defNodes[i]->point());
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
  FvK = Y*Ravg*Ravg/KC;
  cout << " WshellBody          = " << shellBody->totalStrainEnergy() << endl;
  cout << " WpotentialBody_1Ring= " << potentialBody_FirstRing.totalStrainEnergy() << endl;
  cout << " WpotentialBody_2Ring= " << potentialBody_SecondRing.totalStrainEnergy() << endl;
  cout << " Total Energy        = " << (shellBody->totalStrainEnergy() + potentialBody_FirstRing.totalStrainEnergy() + potentialBody_SecondRing.totalStrainEnergy()) << endl;
  cout << " FvK number          = " << FvK              << endl;
  cout << " Asphericity         = " << asphericity      << endl;
  cout << " Ravg                = " << Ravg             << endl;
  cout << "-----------------------------------------"   << endl;
  
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
  





  
  // Clean up
  delete shellBody;
  
  for (uint i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }
  
  return (0);  
}
