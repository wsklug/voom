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
#include "SCElastic.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "LoopShell.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"

#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "VoomMath.h"

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
  
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  double mu = 1.0;
  double kS = 1.0;
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
  inp >> temp >> KC;
  inp >> temp >> KG;
  inp >> temp >> C0;
  inp >> temp >> mu;
  inp >> temp >> kS;
  inp >> temp >> quadOrder;


  
  inp.close();

  // List input parameters
  cout << " modelName               : " << modelName         << endl
       << " KC                      : " << KC                << endl
       << " KG                      : " << KG                << endl
       << " C0                      : " << C0                << endl
       << " mu                      : " << mu                << endl
       << " kS                      : " << kS                << endl
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
  vector<DeformationNode<3>* > defNodes;
   
  
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
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    
    NodeBase::DofIndexMap idx(3);
    for(uint j = 0; j < 3; j++) idx[j] = dof++;
    
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }

   

  // Read in triangle connectivities
  vector<tvmet::Vector<int,3> > connectivities;
  vector<vector<int > > memconn;
  tvmet::Vector<int,3> ct;
  vector<int > cm(3,0);
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
    cm[0] = ct(0); cm[1] = ct(1); cm[2] = ct(2);
    // Compute initial average edge length - used in potential body
    AverageEdgeLength += tvmet::norm2(defNodes[ct(0)]->point()- defNodes[ct(1)]->point());
    AverageEdgeLength += tvmet::norm2(defNodes[ct(1)]->point()- defNodes[ct(2)]->point());
    AverageEdgeLength += tvmet::norm2(defNodes[ct(2)]->point()- defNodes[ct(0)]->point());
    
    connectivities.push_back(ct);
    memconn.push_back(cm);
  }
 
  AverageEdgeLength /= double(ntri*3);
  
  cout << "Number of nodes     = " << nodes.size()          << endl
       << "Number of elements  = " << connectivities.size() << endl
       << "Average Edge Length = " << AverageEdgeLength     << endl;

  // Close mesh file in vtk format
  ifs.close();
  


  // -------------------------------------------------------------------
  // Create model
  // -------------------------------------------------------------------
  // Create Loop shell body with bending energy
  SCElastic bendingMaterial(KC,KG,C0);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  LoopShellBody<SCElastic> * bendingBody = new LoopShellBody<SCElastic>(bendingMaterial, 
									 connectivities, 
									 nodes, 
									 quadOrder);
  // Initialize body
  bendingBody->SetRefConfiguration(AverageEdgeLength);
  bendingBody->setOutput(paraview);
  bendingBody->compute(true, false, false);
  cout << "Initial shell body energy = " << bendingBody->totalStrainEnergy() << endl;
  cout << "Initial shell body energy = " << bendingBody->energy() << endl;
  cout << "Average edge length from body = " << bendingBody->AverageEdgeLength() << endl << endl;
  // shellBody->checkConsistency();
  

  // Create body with in-plane elastic energy
  EvansElastic stretchMaterial(0.0, 0.0, 0.0, mu, kS);

  // Use generic C0MembraneBody for stretching with the EvansElastic material
  C0MembraneBody<TriangleQuadrature, EvansElastic, ShapeTri3> * membraneBody =
    new C0MembraneBody<TriangleQuadrature, EvansElastic, ShapeTri3>(stretchMaterial, 
  								    memconn, 
								    nodes, 
  								    quadOrder);
  membraneBody->setOutput(paraview);
  membraneBody->compute(true, false, false);
  cout << "Initial membrane body energy = " << membraneBody->totalStrainEnergy() << endl;
  cout << "Initial membrane body energy = " << membraneBody->energy() << endl;
  // membraneBody->checkConsistency();
  



  // Create Model and solver
  Model::BodyContainer bdc;
  bdc.push_back(bendingBody);  
  bdc.push_back(membraneBody); 
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

  // Print before everything
  char vtkname[100]; 
  sprintf(vtkname, "Capsid-0");
  bendingBody->printParaview(vtkname);

  // Solve
  solver.solve(&model);

  // Print after solving
  sprintf(vtkname,"Capsid-1");
  bendingBody->printParaview(vtkname);


  // Hard coded print of outline
  ofstream ofs1("T4outSolved.vtk");
  if (!ofs1) {
    cout << "Error: cannot open outline file." << endl;
    exit(0);
  }
  
  ofs1 << "# vtk DataFile Version 3.0" << endl
	 << "vtk output" << endl
	 << "ASCII" << endl
 	 << "DATASET POLYDATA" << endl
         << "POINTS " << 182 << " float" << endl;

  for (uint i = 0; i < 182; i++)
  {
    const Vector3D & nodalPos =  defNodes[i]->point();
    ofs1 << std::setprecision(16) 
	 << nodalPos(0) << "  "
	 << nodalPos(1) << "  "
	 << nodalPos(2) << endl;
  } 

  // Outline for capsomers
  ifstream ifs2("T4out.vtk");
  if (!ifs) {
    cout << "Error: Cannot open input file T4out.vtk " << endl;
    exit(0);
  }

  ifs2 >> token; 
  while( token != "LINES") ifs2 >> token; 
  ofs1 << token << " ";
  ifs2 >> token;
  ofs1 << token << " ";
  ifs2 >> token;
  ofs1 << token << endl;
  ifs2 >> token;
  while(!ifs2.eof())
  {
    ofs1 << token << " ";
    ifs2 >> token;
  }
  ifs2.close();
  // End of hard coded print of outline




  bendingBody->compute(true, false, false);
  membraneBody->compute(true, false, false);
  cout << "Energy in bending body  = " << bendingBody->energy() << endl;
  cout << "Energy in membrane body = " << membraneBody->energy() << endl;
  cout << "Total energy            = " << bendingBody->energy() + membraneBody->energy() << endl;
  
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
  

  



  // -------------------------------------------------------------------
  // Clean up
  // -------------------------------------------------------------------
  delete bendingBody;
  delete membraneBody;
  
  for (uint i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }
  
  return (0);  
}
