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
  double SearchRadius = 0.0;
  double MaxForce = 1.0;
  double nForceSteps = 1.0;
  // Potential input parameters


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
  inp >> temp >> SearchRadius;
  inp >> temp >> MaxForce;
  inp >> temp >> nForceSteps;

  
  inp.close();

  // List input parameters
  cout << " modelName     : " << modelName    << endl
       << " OutlineName   : " << OutlineName  << endl
       << " min           : " << min          << endl
       << " KC            : " << KC           << endl
       << " KG            : " << KG           << endl
       << " C0            : " << C0           << endl
       << " mu            : " << mu           << endl
       << " kS            : " << kS           << endl
       << " NPout         : " << NPout        << endl
       << " OutFile       : " << OutFile      << endl
       << " Refinement    : " << refinement   << endl
       << " Search Radius : " << SearchRadius << endl
       << " Max Force     : " << MaxForce     << endl
       << " N force steps : " << nForceSteps  << endl;
      
  





  // Read in the mesh, in (ascii) legacy vtk format.
  // Create input stream
  ifstream ifs;
  ifs.open(modelName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << modelName << endl;
    exit(0);
  }
 
  // Create vector of nodes
  unsigned int dof = 0, npts = 0;
  vector<NodeBase* > nodesShell, nodesSolver;
  vector<DeformationNode<3>* > defNodes, ForceNodesV0, ForceNodesV1;
  double Ravg = 0;

  // Opposite vertices where forces are applied
  double phi = 0.5*(1 + sqrt(5.0));
  Vector3D V0(0.0, 1.0,-phi);
  Vector3D V1(0.0,-1.0, phi);
  uint CountV0nodes = 0, CountV1nodes = 0, CountSolverNodes = 0;
   
  
  // Input .vtk file containing nodes and connectivities
  string token;
  while( token != "POINTS" ) ifs >> token;
  ifs >> npts; 
  nodesSolver.reserve(npts);
  nodesShell.reserve(npts);
  defNodes.reserve(npts);
  ifs >> token;   // skip number type

  // read in points
  // temporary nodal position container
  vector<DeformationNode<3>::Point > xlist, xlistV0, xlistV1;
  for(uint i = 0; i < npts; i++) {
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    
    if (tvmet::norm2(x-V0) < SearchRadius) {
      xlistV0.push_back( x ); // Nodes where force is applied in V10 direction
    }
    else if (tvmet::norm2(x-V1) < SearchRadius) {
      xlistV1.push_back( x ); // Nodes where force is applied in V01 direction
    }
    else {
      xlist.push_back( x ); // free nodes
    }
  } // end of reading points
  CountV0nodes = xlistV0.size();
  CountV1nodes = xlistV1.size();
  CountSolverNodes = xlist.size();

  // Now create solver nodes
  for(uint i = 0; i < CountSolverNodes; i++) {
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) idx[j] = dof++;
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,xlist[i]);
    nodesShell.push_back( n );
    nodesSolver.push_back( n )
    defNodes.push_back( n );
  }
  for(uint i = 0; i < CountV0Nodes; i++) {
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) idx[j] = dof++;
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,xlist[i]);
    nodesShell.push_back( n );
    nodesSolver.push_back( n )
    defNodes.push_back( n );
  }
  
  Ravg /= nodesShell.size();
  cout << endl << "V0 nodes number = " << CountV0nodes << endl;
  cout << endl << "V1 nodes number = " << CountV1nodes << endl;
  assert(CountV0nodes == CountV1nodes);
  
   
  // Read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  int ntri = 0, tmp = 0;
  ifs >> ntri;
  connectivities.reserve(ntri);
  
  ifs >> tmp;
  for (uint i = 0; i < ntri; i++)
  {
    ifs >> tmp;
    if(tmp != 3) cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    for(uint a = 0; a < 3; a++) ifs >> ct[a];
    connectivities.push_back(ct);
  }
  cout << "Number of nodes    = " << nodesShell.size()          << endl
       << "Ravg               = " << Ravg                  << endl;
  cout << "Number of elements = " << connectivities.size() << endl;
  
  // Close mesh file in vtk format
  ifs.close();
  
 



  

  // Create Loop shell body with both bending and in-plane elastic energy
  int quadOrder = 1;
  
  // Bending material object, to be copied when body generates new elements
  EvansElastic bodyMaterial(KC, KG, C0, mu, kS);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  LoopShellBody<EvansElastic> * shellBody = new LoopShellBody<EvansElastic>(bodyMaterial, connectivities, nodesShell, quadOrder);
  // Initialize body
  shellBody->compute(true, false, false);
  cout << "Initial shell body energy = " << shellBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << shellBody->energy() << endl;
        // bending_body->checkConsistency();

  




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
  
  Model model(bdc, nodesSolver);
       // Consistency check
       // model.checkConsistency(true,false);
       // model.checkRank(model.dof()-6,true);
  
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
  // Set bounds on eta variables (if any)
  // if ((min == 1 || min ==3) && WcType == 6)
  // {
  //   solver.setBounds(BoundType, LowerBound, UpperBound);
  // }
  // cout << "Bound  DOF = " << NumDOF << endl;
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

   // Apply force
   double ForceStep = MaxForce/nForceSteps/CountV1nodes;
   Vector3D V01(0.0), Force(0.0);
   V01 = V1-V0; V01 = V01/tvmet::norm2(V01);
   for (uint step = 0; step < nForceSteps; step++) {
     // Apply force
     for (uint count = 0; count < ForceNodesV1.size(); count++) {
       Force = V01*ForceStep;
       ForceNodesV1[count]->updateForce(Force); // update force add the force step to the nodal forces
     }
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
   cout << "WshellBody           = " << shellBody->totalStrainEnergy() << endl;
   cout << " FvK number          = " << FvK              << endl;
   cout << " Asphericity         = " << asphericity      << endl;
   cout << " Ravg                = " << Ravg             << endl;
   cout << "-----------------------------------------"   << endl;
  
   time (&end);
   dif = difftime (end,start);
   cout << endl << "All done :) in " << dif  << " s" << endl;
 





  
   // Clean up
   delete shellBody;

   for (uint i = 0; i<nodesShell.size(); i++)
   {
     delete nodesShell[i];
   }

   return (0);  
}
