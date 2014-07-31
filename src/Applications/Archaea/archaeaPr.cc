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
#include "ProteinBody.h"
#include "ProteinLennardJones.h"
#include "MontecarloProtein.h"
#include "LoopShell.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "VoomMath.h"
#include "Utils/PrintingProtein.h"

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
  string initialConf;
  int ICprovided = 0;
  double ICtol = 1.0e-5;
  string outputFileName;
  int VTKflag = 0;
  double InflationFactor = 1.0;
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  int quadOrder = 1;
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

  // Monte Carlo solver parameters
  int nMCsteps = 1000;
  double T01, T02, FinalRatio;
  int CompNeighInterval;
  int MCmethod = 0;

  // Number of nodes per protein
  int NnodePr = 10;
  int PrMaxNum = 10;
  int PrintEvery = 10;


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
  inp >> temp >> initialConf;
  inp >> temp >> ICprovided;
  inp >> temp >> ICtol;
  inp >> temp >> VTKflag;
  inp >> temp >> InflationFactor;
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> quadOrder;
  inp >> temp >> penaltyVolumeInit;
  inp >> temp >> penaltyAreaInit;
  inp >> temp >> penaltyVolumeTol;
  inp >> temp >> penaltyAreaTol;
  inp >> temp >> PotentialSearchRF;
  inp >> temp >> epsilon;
  inp >> temp >> sigma;
  inp >> temp >> Rshift;
  inp >> temp >> RconnSF;
  inp >> temp >> nMCsteps;
  inp >> temp >> T01;
  inp >> temp >> T02;
  inp >> temp >> FinalRatio;
  inp >> temp >> CompNeighInterval;
  inp >> temp >> MCmethod;
  inp >> temp >> NnodePr;
  inp >> temp >> PrMaxNum;
  inp >> temp >> PrintEvery;


  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " outputFileName          : " << outputFileName    << endl
       << " initialConfiguration    : " << initialConf       << endl
       << " IC provided             : " << ICprovided        << endl
       << " IC tol                  : " << ICtol             << endl
       << " VTKflag                 : " << VTKflag           << endl
       << " InflationFactor         : " << InflationFactor   << endl
       << " KC                      : " << KC                << endl
       << " KG                      : " << KG                << endl
       << " C0                      : " << C0                << endl
       << " Loop Shell quadOrder    : " << quadOrder         << endl
       << " Penalty volume factor   : " << penaltyVolumeInit << endl
       << " Penalty area factor     : " << penaltyAreaInit   << endl  
       << " Penalty volume tol      : " << penaltyVolumeTol  << endl
       << " Penalty area tol        : " << penaltyAreaTol    << endl  
       << " Potential Search factor : " << PotentialSearchRF << endl
       << " epsilon                 : " << epsilon           << endl
       << " sigma                   : " << sigma             << endl
       << " Rshift                  : " << Rshift            << endl
       << " Rconn search radius     : " << RconnSF           << endl
       << " MC steps                : " << nMCsteps          << endl
       << " MC T01                  : " << T01               << endl
       << " MC T02                  : " << T02               << endl
       << " MC FinalRatio           : " << FinalRatio        << endl
       << " MC ComputeNeighInterval : " << CompNeighInterval << endl
       << " MC method               : " << MCmethod          << endl
       << " Number of nodes per Pr  : " << NnodePr           << endl
       << " Max number of Pr        : " << PrMaxNum          << endl
       << " Print .vtk file every   : " << PrintEvery        << endl;



      
  






  // -------------------------------------------------------------------
  // Read mesh and setup geometric model
  // -------------------------------------------------------------------
  // Read in the mesh and, if provided the initial configuration
  vector<DeformationNode<3>::Point > IC;
  unsigned NumIC = 0;
  vector<uint > ICcheck;
  if (ICprovided == 1) {
    ifstream ifsIC;
    ifsIC.open(initialConf.c_str(), ios::in);
    if (!ifsIC) {
      cout << "Cannot open input file: " << initialConf << endl;
      exit(0);
    }
    
    string line;
    ifsIC >> line;
    ifsIC >> NumIC;
    ICcheck.assign(NumIC, 1);
    for(uint i = 0; i < NumIC; i++) {
      uint temp = 0;
      DeformationNode<3>::Point xIC;
      if (VTKflag == 1) {
	ifsIC >> xIC(0) >> xIC(1) >> xIC(2); 
      }
      else {
	ifsIC >> temp >> xIC(0) >> xIC(1) >> xIC(2); }
      IC.push_back(xIC);
    }
  } // end of ICprovided loop
  cout << "NumIC = " << IC.size() << endl;

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
  vector<DeformationNode<3>* > defNodes;

  // Input .vtk or inp files containing nodes and connectivities
  string token;
  ifs >> token;
  ifs >> npts; 
  NumDoF = npts*3; // Assumed all nodes have 3 DoF
  nodes.reserve(npts);
  defNodes.reserve(npts);
  vector<ProteinNode *> Proteins;

  // read in points
  uint NextProtein = 0;
  for(uint i = 0; i < npts; i++) {
    uint NodeNum = 0;
    DeformationNode<3>::Point x;
    if (VTKflag == 1) {
      ifs >> x(0) >> x(1) >> x(2); }
    else {
      ifs >> NodeNum >> x(0) >> x(1) >> x(2); }
    
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) { 
      idx[j] = dof++;
    }
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
    
    if (ICprovided == 1) {
      if (VTKflag == 1) {
	if (i < NumIC) {
	  ProteinNode * PrNode = new ProteinNode(n);
	  Proteins.push_back(PrNode);
	}
      }
      else {
	for (uint k = 0; k < NumIC; k++) {
	  if (tvmet::norm2(IC[k] - x) < ICtol && ICcheck[k] == 1) {
	    ProteinNode * PrNode = new ProteinNode(n);
	    Proteins.push_back(PrNode);
	    ICcheck[k] = 0;
	    break;
	  }
	}
      }
    }
    else {
      if (i == NextProtein && Proteins.size() < PrMaxNum) { // create a protein every NodePr nodes and additional conditions
	// if (x(2) > 1.0) { // just for cone example - start from proteins packed close to the tip
	ProteinNode * PrNode = new ProteinNode(n);
	Proteins.push_back(PrNode);
	NextProtein += (rand() % NnodePr) + 1;
	// }
      }
    }

    // Inflate vessel if necessary
    x *= InflationFactor;
    n->setPoint(x);
    n->resetPosition();
      
    nodes.push_back( n );
    defNodes.push_back( n );
   
  } // end of reading points

  ICcheck.clear();
   
  // Read in triangle connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  uint ntri = 0, ElemNum = 0, tmp = 0;
  map<DeformationNode<3> *, vector<DeformationNode<3> *> > PossibleHosts;

  ifs >> token;
  ifs >> ntri;
  connectivities.reserve(ntri);
  for (uint i = 0; i < ntri; i++)
  {
    ifs >> ElemNum;
    ifs >> ct(0) >> ct(1) >> ct(2);
    if (VTKflag == 0) {
      ct(0) -= 1;  ct(1) -= 1;  ct(2) -= 1; }
      
    connectivities.push_back(ct);
  } // end of reading triangle connectivities

  // build possible hosts table
  for(uint i = 0; i < npts; i++)
  {
    set<DeformationNode<3> *> setHosts;
    for (uint j = 0; j < ntri; j++) {
      ct = connectivities[j];
      if (i == ct(0) || i == ct(1) || i == ct(2) ) {
	setHosts.insert( defNodes[ct(0)] );
	setHosts.insert( defNodes[ct(1)] );
	setHosts.insert( defNodes[ct(2)] );
      }
    }
    
    vector<DeformationNode<3> *> vectorHosts;
    for (set<DeformationNode<3> *>::iterator pHost = setHosts.begin();
	 pHost != setHosts.end(); pHost++) {
      vectorHosts.push_back(*pHost);
    }

    PossibleHosts[ defNodes[i] ] = vectorHosts;
  }
  
  cout << "Number of nodes     = " << nodes.size()          << endl
       << "Number of proteins  = " << Proteins.size()       << endl
       << "Number of elements  = " << connectivities.size() << endl;

  // Close mesh file in inp format
  ifs.close();
  


  // -------------------------------------------------------------------
  // Create model
  // -------------------------------------------------------------------
  // Create Loop shell body with both bending and in-plane elastic energy
  /*
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
  // shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
  cout << "Initial shell body volume = " << shellBody->volume() << endl;
  cout << "Initial shell body area   = " << shellBody->area()   << endl;
  cout << "Average edge length from body = " << shellBody->AverageEdgeLength() << endl << endl;
  // shellBody->checkConsistency();
  
  */

  // Initiliaze potential material
  ProteinLennardJones Mat(epsilon, sigma);

  // Then initialize potential body
  ProteinBody * PrBody = new ProteinBody(Proteins, &Mat, PotentialSearchRF);
  
  PrBody->compute(true, false, false);
  cout << "Initial protein body energy = " << PrBody->energy() << endl;






  /*
  // Create Model and solver
  Model::BodyContainer bdc;
  bdc.push_back(shellBody);
  Model model(bdc, nodes);
       // Consistency check
       // model.checkConsistency(true,false);
       // model.checkRank(model.dof()-3,true);
  */
  
  // Initialize printing utils
  cout << "Seraching connectivity over R = " <<  RconnSF << endl;
  PrintingProtein PrintArchaea(modelName, outputFileName+"iter", nodes, connectivities, Proteins, RconnSF);
  
  // Initialize Montecarlo Solver
  MontecarloProtein MCsolver(Proteins, PrBody, PossibleHosts, MCmethod, &PrintArchaea, PrintEvery, nMCsteps);
  MCsolver.SetTempSchedule(MCsolver.EXPONENTIAL, T01, T02, FinalRatio);
  MCsolver.solve(CompNeighInterval, PotentialSearchRF);

  /*
  // Print proteins positions for further processing
  for (uint i = 0; i < Proteins.size(); i++)
  {
    DeformationNode<3>::Point P =  Proteins[i]->getHostPosition();
    cout << P[0] << " " << P[1] << " " << P[2] << endl;
  }
  */


  // End of program
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
  

  
  // -------------------------------------------------------------------
  // Clean up
  // -------------------------------------------------------------------
  // delete shellBody;
  delete PrBody;
  
  for (uint i = 0; i<Proteins.size(); i++)
  {
    delete Proteins[i];
  }

  for (uint i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }
  
  return (0);  
}
