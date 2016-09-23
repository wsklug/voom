#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#include "voom.h"
#include "NodeBase.h"
#include "Node.h"
#include "LMEtetBody.h"

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
  
  // LME parameters
  double beta = 0.0;           // LME beta
  double search_radius = 0.0;  // LME search radius
  unsigned int quad_order = 1;
  double tol = 1.0e-12;
  unsigned int nItMax = 20;
   
  // Eigenvectors flag
  unsigned int ModesNumber = 10;
  unsigned int EgvFlag = 0;

  // Material parameters
  double rho = 1.0;
  // rho=17298.712055/bd.volume()*1.66e-27;   //density calculated from atomic mass
  double E = 200.0;
  double nu = 0.4;
   
  // Reading input from file passed as argument
  ifstream inpFile;
  inpFile.open(parameterFileName.c_str(), ios::in);
  if (!inpFile) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  string dump;
  
  inpFile >> dump >> modelName;
  inpFile >> dump >> beta;
  inpFile >> dump >> search_radius;
  inpFile >> dump >> quad_order;
  inpFile >> dump >> tol;
  inpFile >> dump >> nItMax;
  inpFile >> dump >> ModesNumber;
  inpFile >> dump >> EgvFlag;
  inpFile >> dump >> rho;
  inpFile >> dump >> E; 
  inpFile >> dump >> nu;
  
  inpFile.close(); 
 
  // List input parameters
  cout << " modelName        : " << modelName     << endl
       << " beta             : " << beta          << endl
       << " Search radius    : " << search_radius << endl
       << " Quadrature order : " << quad_order    << endl
       << " LME tolerance    : " << tol           << endl
       << " LME max iter     : " << nItMax        << endl
       << " ModesNumber      : " << ModesNumber   << endl
       << " Eigenvector flag : " << EgvFlag       << endl
       << " rho              : " << rho           << endl
       << " E                : " << E             << endl
       << " nu               : " << nu            << endl;


 


  // Read nodal coordinates and tet connectivity
  ifstream LMEmeshFile;
 
  // Create input stream
  LMEmeshFile.open( modelName.c_str(), ios::in);
  if (!LMEmeshFile ) {
    cout << "Cannot open mesh file: " << LMEmeshFile << endl;
    exit(0);
  }

  // Create vector of nodes and connectivity table
  int dof = 0, npts = 0, i = 0, j = 0;
  vector< NodeBase* > nodes;
  vector< DeformationNode<3>* > defNodes;
  vector<vector<unsigned int > > connT;
  
  string token;
  LMEmeshFile >> token; 
  LMEmeshFile >> token; 
  while( token != "POINTS" ) {
    LMEmeshFile >> token;
  }
  LMEmeshFile >> npts; 
  LMEmeshFile >> token; 
  nodes.reserve(npts); 
  defNodes.reserve(npts);

  // read in points
  DeformationNode<3>::Point x;
  for(int i = 0; i < npts; i++)
  {
    LMEmeshFile >> x(0) >> x(1) >> x(2);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  cout << "Number of nodes: " << nodes.size() << endl;

  while( token != "CELLS" ) {
    LMEmeshFile >> token;
  }
  vector<unsigned int > SingleConnT(4,0);
  unsigned int ntri = 0, temp = 0;
  LMEmeshFile >> ntri;
  LMEmeshFile >> token;
  connT.reserve(ntri);
  cout << "Number of elements: " << ntri << endl;

  for (i = 0; i < ntri; i++)
  {
    LMEmeshFile >> temp;
    if(temp != 4) {
      cout << i << " Some mistake reading the elements connectivity from file. Check again." << endl;
    }
    for(j = 0; j < 4; j++) {
      LMEmeshFile >>  SingleConnT[j];
    }
    connT.push_back(SingleConnT);
  }
  LMEmeshFile.close();





  // Create Body
  LMEtetBody protein(defNodes, connT, beta, search_radius, quad_order, tol, nItMax, rho, E, nu);
  std::cout << "Calculating normal modes now." << std::endl;
  unsigned int N = npts*3;
  double* K = new double[N*N];
  for (i = 0; i < N*N; i++) {
    K[i] = 0.0;
  }
  double* eval = new double[N];
  protein.computeNormalModes(K, eval, ModesNumber, EgvFlag);
  for (i = 0; i < N; i++) {
    cout << "i = " << i << " eval[i] = " << eval[i] << endl;
  }




  // protein.checkConsistency(true,false);
  // protein.checkRank(model.dof()-6,true);

  return 0;
}
