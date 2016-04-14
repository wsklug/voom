#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <tvmet/Vector.h>

#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "LoopShellBody.h"
#include "LoopShell.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "VoomMath.h"

using namespace voom;
using namespace std;

void printOutlineParaview(const std::string OutlineName,
			  const int iter,
		          const unsigned int NPout, 
		          const Model::BodyContainer & bdc);

int main(int argc, char* argv[])
{
  time_t start, end;
  double dif;
  time (&start);

  if(argc < 2){
    cout << "Input file missing." << endl;
    return(0);
  }
 
  const double PI = 3.14159265;
  
  // Material properties
  double KC = 1.0;
  double nu = 0.3;
  double C0 = 0.0;
  double Ymin = 1.0;
  double Ymax = 100.0;
  int    Ysteps = 10;
  int    quadOrder = 1;
  // Model names
  string modelName;
  string outlineName;
  unsigned int NPout = 0; // Nodes needed to build the outline

  // Input output
  string parameterFileName = argv[1];
  ifstream inp, ifs;
 
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  inp >> modelName >> KC >> nu >> C0 >> Ymin >> Ymax >> Ysteps >> quadOrder >> outlineName >> NPout;
  inp.close();

  // List input parameters
  cout << " modelName  : " << modelName   << endl
       << " KC         : " << KC          << endl
       << " nu         : " << nu          << endl
       << " C0         : " << C0          << endl
       << " Ymin       : " << Ymin        << endl
       << " Ymax       : " << Ymax        << endl
       << " Ysteps     : " << Ysteps      << endl
       << " quadOrder  : " << quadOrder   << endl
       << " outlineName: " << outlineName << endl
       << " NPout      : " << NPout       << endl;
  


  // Read in the mesh, in (ascii) legacy vtk format.
  string meshFileName = modelName + ".vtk";

  // Create input stream
  ifs.open( meshFileName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << meshFileName << endl;
    exit(0);
  }

  // Create vector of nodes
  int dof = 0, id = 0, npts = 0;
  vector< NodeBase* > nodes;
  vector< DeformationNode<3>* > defNodes;
  double Ravg = 0;
  char key;  ifs >> key;
  
  // Input .vtk file containing nodes and connectivities
  string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  ifs >> npts;    // Read in number of nodes
  cout << "npts = " << npts << endl;
  defNodes.reserve(npts);
  ifs >> token;   // skip number type
  
  // Read in points
  for(int i = 0; i < npts; i++) {
    id = i;
    DeformationNode<3>::Point X;
    ifs >> X(0) >> X(1) >> X(2);
    Ravg += tvmet::norm2(X);
    NodeBase::DofIndexMap idx(3);

    for(int j = 0; j < 3; j++) idx[j] = dof++;

    DeformationNode<3>* n = new DeformationNode<3>(id,idx,X);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size() != 0);
  Ravg /= nodes.size();
  cout << "Number of nodes = " << nodes.size() << endl
       << "Ravg            = " << Ravg         << endl;
  
  // Read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  int ntri = 0, tmp = 0;
  double AvgEdgeLength = 0.0;
  ifs >> ntri; // Read in number of elements
  cout << "Number of triangles: " << ntri << endl;
  connectivities.reserve(ntri);
 
  ifs >> tmp; // cout << " tmp = " << tmp << endl;
  for (int i = 0; i < ntri; i++)
  {
    tmp = 0;
    ifs >> tmp;
    if(tmp != 3) cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    for(int a = 0; a < 3; a++) ifs >> ct[a];
    connectivities.push_back(ct);

    // Edge vectors in current config.
    tvmet::Vector<double,3> e31( defNodes[ct[0]]->point() - defNodes[ct[2]]->point() ), 
                            e32( defNodes[ct[1]]->point() - defNodes[ct[2]]->point() ),
                            e12( defNodes[ct[1]]->point() - defNodes[ct[0]]->point() );

    // Add averate edge length for each triangle
    AvgEdgeLength += (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;
  }
  AvgEdgeLength /= ntri;
  cout << "Average edge length = " << AvgEdgeLength << endl;
   
    

  

  // Create bending and stretching body together

  // Bending
  double KG = (nu-1)*KC;
  
  // Bending material object, to be copied when body generates new elements
  typedef EvansElastic ArchaeaMaterial;

  double DeltaY = (Ymax - Ymin)/(double(Ysteps)-1.0);
  for (int s = 0; s < Ysteps; s++) {
    double Y = Ymax - double(s)*DeltaY;
    // Stretching
    double mu = Y/(2.0*(1.0 + nu));
    double kS = Y/(2.0*(1.0 - nu));

    ArchaeaMaterial aMaterial(KC, KG, C0, mu, kS);

    // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
    typedef LoopShellBody<ArchaeaMaterial > LSB;
    LSB * ArchaeaBody = new LSB(aMaterial, connectivities, nodes, quadOrder);
  
    ArchaeaBody->compute(true, false, false);
    cout << "W initial = " << ArchaeaBody->energy() << endl;
    // cout << endl << "Checking Bending body consistency before imposing ref configuration " << endl;
    // bending_body->checkConsistency(); 
  
    // Set reference configuration
    ArchaeaBody->SetRefConfiguration(AvgEdgeLength);
    ArchaeaBody->compute(true, false, false);
    cout << "W after setting AvgEdgeLength = " << ArchaeaBody->totalStrainEnergy() << endl;
    // cout << endl << "Checking Bending body consistency after imposing ref configuration " << endl;
    // bending_body->checkConsistency();
  
    // Create Model
    Model::BodyContainer bdc;
    bdc.push_back(ArchaeaBody);
    Model model(bdc, nodes);
    // Consistency check
    // model.checkConsistency(true,false);
    // model.checkRank(model.dof()-6,true);

    //  double FvK = 0.0, Y = 0.0;
    // Y = 4.0*kS*mu/(kS+mu);
    // FvK = Y*Ravg*Ravg/KC;
    // cout << "FvK number before solving = " << FvK << endl;

    // Set solver and solve
    int iprint = 0;
    double factr = 1.0e1;
    double pgtol = 1.0e-8;
    int m = 5;
    int maxIter = 1000000;
    std::ifstream lbfgsbinp("lbfgsb.inp");
    lbfgsbinp >> iprint >> factr >> pgtol >> m;
    cout << "Input iprint: " << iprint << endl
	 << "Input factr:  " << factr  << endl
	 << "Input pgtol:  " << pgtol  << endl
	 << "Input m:      " << m      << endl;
  
    Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);
    // Solve to minimize the energy
    // bdc[0]->printParaview("Initial");

    solver.solve(&model);
    
    stringstream OutA;
    OutA << "ArchaeaMin_" << Ysteps - s;
    string OutFileNameA = OutA.str();
    bdc[0]->printParaview(OutFileNameA);
    printOutlineParaview(outlineName, Ysteps - s, NPout, bdc);

    //---------------//
    // Print results //
    cout << "W final = " << ArchaeaBody->energy()                << endl;
    cout << "Area    = " << ArchaeaBody->area()                  << endl;
    cout << "Volume  = " << ArchaeaBody->volume()                << endl;
    cout << "FvKnew  = " << Y*ArchaeaBody->area()/KC             << endl;
    cout << "TotCurv = " << ArchaeaBody->totalCurvature()        << endl;
    cout << "EqAsph  = " << ArchaeaBody->equivalentAsphericity() << endl;

    // Clean up
    delete ArchaeaBody;

  }
  
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
  


  // Clean up
  // delete ArchaeaBody;

  for (int i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }

  return (0);
  
}



void printOutlineParaview(const std::string OutlineName,
			  const int iter,
		          const unsigned int NPout, 
		          const Model::BodyContainer & bdc)
{
  // OutPut bending and stretching energy in current configurationoutline in current configuration

  std::string OutlineFile = OutlineName+".vtk";
  stringstream OutA;
  OutA << OutlineName << "New_" << iter << ".vtk";
  string OutNewFile = OutA.str();
  ofstream ofs(OutNewFile.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
	      << OutNewFile
	      << ") file." << std::endl;
    exit(0);
  }
  unsigned int ind = 0;

  // Extract elements
  Body::ElementContainer Elements = bdc[0]->elements();
  Body::NodeContainer Nodes = bdc[0]->nodes();

  // Node Section
  ofs << "# vtk DataFile Version 3.0" << endl
      << "vtk output" << endl
      << "ASCII" << endl
      << "DATASET POLYDATA" << endl
      << "POINTS " << NPout << " float" << endl;

  // Output nodal postions
  for (int n = 0; n < NPout; n++)
  {
    DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(Nodes[n]);

    if (node == NULL) exit(0);

    const Vector3D & nodalPos =  node->point();
    ofs << std::setprecision(16) 
	<< nodalPos(0) << "  "
	<< nodalPos(1) << "  "
	<< nodalPos(2) << std::endl;
  }

  //   // Outline for pentamers, hexamers, and heptamers
  ifstream ifs(OutlineFile.c_str());
  if (!ifs) {
    cout << "Cannot open input file: " << OutlineFile << endl;
    exit(0);
  }

  string token;
  ifs >> token; 
  while( token != "LINES") ifs >> token; 
  ofs << token << " ";
  ifs >> token;
  ofs << token << " ";
  ifs >> token;
  ofs << token << endl;
  ifs >> token;
  while(!ifs.eof())
  {
    ofs << token << " ";
    ifs >> token;
  }
  ifs.close();

  ofs.close();

}





    
    
   
