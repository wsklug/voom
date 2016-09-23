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
#include "C0MembraneBody.h"
#include "StretchHexonBody.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "MontecarloTwoStages.h"
#include "VoomMath.h"

#include "Utils/PrintingStretches.h"

// #include "print-body.cc"

using namespace voom;
using namespace std;

/*! This program creates a shell model with LoopShellBody for bending
  and 2 C0MembraneBodies for stretching, one of pentamers and one of
  hexamers.
 */

int main(int argc, char* argv[])
{
  time_t start,end;
  double dif;
  time (&start);

  int i = 0, j = 0;
  if(argc < 2){
    cout << "Input file missing." << endl;
    return(0);
  }
 
  const double PI = 3.14159265;
  
  int min = 1;           // minimize with respect to eta flag (1-> do it, 0-> don't do it)
  int Ns = 1;            // Number of shear dof
  double in_eta = 0.0;   // eta initial guess
  double in_theta = 0.0; // theta initial guess
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  // Stretching
  double mu = 1.0;
  double kS = 1.0;
  double AngleShift = 0.0;
  double HexonShift = 0.0;
  int WcType = 1;
  double WcConst0 = 1.0;
  double WcConst1 = 1.0;
  double WcConst2 = 1.105;
  double WcConst[3];

  // Parameters for MonteCarlo solver
  unsigned int maxIterMTS = 100;
  double finalRatioMTS = 1.0e-5;
  double T1 = 300.0;
  double T2 = 1.0;
  int VarEta = -2;
  int VarTheta = 0;
  
  string modelName;
  string outlineName;
  unsigned int NPout = 0; // Nodes needed to build the outline
  string parameterFileName = argv[1];
  ifstream inp, ifs;
 
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  inp >> modelName >> outlineName >> NPout >> in_eta >> in_theta >> KC >> KG >> C0 >> mu >> kS >> min >> Ns >> AngleShift >> HexonShift >> WcType >> WcConst0 >> WcConst1 >> WcConst2 >> maxIterMTS >> finalRatioMTS >> T1 >> T2 >> VarEta >> VarTheta;
  inp.close();

  // List input parameters
  cout << " modelName   : " << modelName     << endl
       << " outlineName : " << outlineName   << endl
       << " NPout       : " << NPout         << endl
       << " in_eta      : " << in_eta        << endl
       << " in_theta    : " << in_theta      << endl
       << " KC          : " << KC            << endl
       << " KG          : " << KG            << endl
       << " C0          : " << C0            << endl
       << " mu          : " << mu            << endl
       << " kS          : " << kS            << endl
       << " min         : " << min           << endl
       << " NstretchDOF : " << Ns            << endl
       << " AngleShift  : " << AngleShift    << endl
       << " HexonShift  : " << HexonShift    << endl
       << " WcType      : " << WcType        << endl
       << " keta        : " << WcConst0      << endl
       << " ktheta      : " << WcConst1      << endl
       << " eta bar     : " << WcConst2      << endl
       << " MaxIter     : " << maxIterMTS    << endl
       << " FinalRatio  : " << finalRatioMTS << endl
       << " T1          : " << T1            << endl
       << " T2          : " << T2            << endl
       << " VarEta      : " << VarEta        << endl
       << " VarTheta    : " << VarTheta      << endl; 
  
  // Populate Wconst
  WcConst[0] = WcConst0;
  WcConst[1] = WcConst1;
  WcConst[2] = WcConst2;
  
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
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;   // skip number type
  cout << "npts = " << npts << endl;
  // read in points
  for(i = 0; i < npts; i++) {
    id = i;
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);

    for(j = 0; j < 3; j++) idx[j] = dof++;

    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  assert(nodes.size() != 0);
  Ravg /= nodes.size();
  cout << "Number of nodes = " <<nodes.size() << endl
       << "Ravg            = " << Ravg        << endl;
  
  // Read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  int ntri = 0, tmp = 0, a = 0;
  ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "Number of triangles: " << ntri << endl;

  ifs >> tmp; // cout << " tmp = " << tmp << endl;
  for (i = 0; i < ntri; i++)
  {
    tmp = 0;
    ifs >> tmp;
    if(tmp != 3) cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    for(a = 0; a < 3; a++) ifs >> ct[a];
    connectivities.push_back(ct);
  }



  // Read in the stretch directions. 
  while( token != "shear_direction") ifs >> token; 
  double AvgEdgeLength = 0.0, stretch = 0.0, TOL = 2.5e-1, stretch_dir_norm = 0.0, factor = 0.0; 
  vector<double> stretch_angle;
  vector< std::vector<int> > hex_connectivities;
  std::vector<int> cm(3,0.0);
  tvmet::Vector<double,3> stretch_dir;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >>  stretch_dir(0) >> stretch_dir(1) >> stretch_dir(2);
    stretch_dir_norm = tvmet::norm2(stretch_dir); 
    
    for(j = 0; j < 3; j++) {
      cm[j] = connectivities[i](j);
    }
    
    // Edge vectors in current config.
    tvmet::Vector<double,3> e31(defNodes[cm[0]]->point()-defNodes[cm[2]]->point()), 
                            e32(defNodes[cm[1]]->point()-defNodes[cm[2]]->point()),
                            e12(defNodes[cm[1]]->point()-defNodes[cm[0]]->point());
    // Compute averate edge length for each triangle
    AvgEdgeLength += (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;

    stretch_dir = stretch_dir/stretch_dir_norm;
    hex_connectivities.push_back(cm);
    // Calculate the stretch angle
      
    //normalize the edge vectors
    e31=e31/tvmet::norm2(e31);
    e32=e32/tvmet::norm2(e32);
    e12=e12/tvmet::norm2(e12);

      // Which edge is closest to the stretch vector?  
      // Take cross product.  If zero (or < TOL) then two vectors are parallel or anti-parallel.  
      // Compute angles in degrees.
    if(tvmet::norm2(tvmet::cross(stretch_dir,e31))<TOL) {

        if(tvmet::dot(stretch_dir,e31)>0) stretch=0.; // parallel with edge 31
        else stretch=180.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(stretch_dir,e32))<TOL){

        if(tvmet::dot(stretch_dir,e32)>0) stretch=60.; // parallel with 32
        else stretch=240.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(stretch_dir,e12))<TOL) {

        if(tvmet::dot(stretch_dir,e12)>0) stretch=120.; // parallel with 12
        else stretch=300.; // anti-parallel

      } 
      // If TOL is too small, it bombs.
      else {std::cout<<"Problem. The stretch direction is not parallel to any triangle edge"<<std::endl; exit(1);}
      
    // Convert from degrees to radians.
    stretch_angle.push_back((stretch+AngleShift)*PI/180.); 
  }
  AvgEdgeLength /= ntri;
  cout << "AvgEdgeLength = " << AvgEdgeLength << endl;
  cout << "Number of hexon elements = "  <<hex_connectivities.size() << endl;

  // Read capsomers to assign different stretch dof to each one
  vector<unsigned int > CapsomersNum(ntri, 0.0);
  while( token != "capsomer") ifs >> token;
  ifs >> token;
  ifs >> token;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >> CapsomersNum[i];
  }
  cout << "Number of capsomers = " <<  CapsomersNum[ntri-1]+1 << endl;
  // Close mesh file in vtk format
  ifs.close();
  


  // Initialize stretch node (ScalarFieldNode type)
  vector<ScalarFieldNode<3>* > stretchNodes, directionNodes;
  stretchNodes.reserve(Ns);   directionNodes.reserve(Ns);
  
  NodeBase::DofIndexMap idE(1);
  ScalarFieldNode<3>::PositionVector pE(0.0);
  for (i = 0; i < Ns; i++)
  {
    id++;
    idE[0] = dof++;
    ScalarFieldNode<3>* stretchN = new ScalarFieldNode<3>(id, idE, pE, in_eta);
    stretchNodes.push_back(stretchN);
  }

  // ! To keep some of the dof constant we can decide to not insert them in Model BUT
  // THIS ONLY WORKS IF the constant dof are numbered last
  for (i = 0; i < Ns; i++)
  {
    id++;
    idE[0] = dof++;
    ScalarFieldNode<3>* dirN = new ScalarFieldNode<3>(id, idE, pE, in_theta);
    directionNodes.push_back(dirN);
  }
  


  // Create bending body; SCElastic only has curvature terms, no in-plane strain terms
  int quadOrder = 1;
  // Bending material object, to be copied when body generates new elements
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);
  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  typedef LoopShellBody<BendingMaterial> LSB;
  LSB bending_body = LSB(bending, connectivities, nodes, quadOrder);





  
  // Hexamers
  // Use specialezed SkewHexonBody (modified version of C0MembraneBody), which has material hard-coded as
  // ModifiedEvansElasticMin, which takes care of the multiplicative F=AG "Eigen-strain" decomposition.
  TriangleQuadrature  HexQuad(1);
  const ShapeTri3::CoordinateArray dummy(0.0);
  ShapeTri3 HexShape(dummy);
 StretchHexonBody * hex_body = new StretchHexonBody(hex_connectivities, defNodes, stretchNodes, directionNodes, stretch_angle, mu, kS, WcType, WcConst, HexQuad, &HexShape, CapsomersNum);
  hex_body->setOutput(paraview);
  hex_body->compute(true, false, false);
  cout << "En_Hex = " << hex_body->totalStrainEnergy() << endl;
  cout << "En_Hex = " << hex_body->energy() << endl;
  cout << "En_Hex = " << hex_body->totalConformationalEnergy() << endl;
  // hex_body->checkConsistency();
  

  
  // Set reference configuration
  hex_body->SetRefConfiguration(AvgEdgeLength);
  cout << "En_Hex = " << hex_body->totalStrainEnergy() << endl;
  cout << "En_Hex = " << hex_body->energy() << endl;
  


  
  


  // Create Model
  Model::BodyContainer bdc;
  bdc.push_back(hex_body);
  if (min == 1)   // shear magnitude is a dof, otherwise energy is minimized keeping eta fixed
  {
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
  }
  else if (min == 2)   // shear magnitude and direction are dof, otherwise energy is minimized keeping eta and theta fixed
  {
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
  }
  Model model(bdc, nodes);
  // Consistency check
  // model.checkConsistency(true,false);
  // model.checkRank(model.dof()-6,true);



  


  
  // Set solver and solve
  int iprint = 0;
  double factr = 1.0e1;
  double pgtol = 1.0e-6;
  int m = 5;
  int maxIter = 10000;
  // std::ifstream lbfgsbinp("lbfgsb.inp");
  // lbfgsbinp >> iprint >> factr >> pgtol >> m;
  cout << "Input iprint: " << iprint << endl
       << "Input factr:  " << factr  << endl
       << "Input pgtol:  " << pgtol  << endl
       << "Input m:      " << m      << endl;
 
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);


  PrintingStretches PrintHexagonalLattice(meshFileName, outlineName, NPout, 
					  "MTiteration", "MToutline", "MTdirections", 0, 
					  bdc, stretchNodes, directionNodes, CapsomersNum, 
					  AngleShift, HexonShift, false);
  cout << "here3" << endl;
  PrintHexagonalLattice.printMaster(-2);
 cout << "here4" << endl;
  PrintHexagonalLattice.setPrintingConfig(true);

  vector<int > VarType;
  vector<vector<ScalarFieldNode<3>* > > MontecarloDof;
  if (VarTheta > -2)
  {
    MontecarloDof.push_back(directionNodes);
    VarType.push_back(VarTheta);
  }
  if (VarEta > -2)
  {
    MontecarloDof.push_back(stretchNodes);
    VarType.push_back(VarEta);
  }
  
  MontecarloTwoStages MTS(MontecarloDof, VarType, &solver, &PrintHexagonalLattice, maxIterMTS, true); 
  MTS.SetTempSchedule(MTS.EXPONENTIAL, T1, T2, finalRatioMTS);
  MTS.solve(&model); 
  
  PrintHexagonalLattice.printMaster(-1);
 
  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
 

  
  // Clean up
  delete hex_body;

  for (i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }

  if (min==0)
  {
    for (i = 0; i<Ns; i++)
    {
      delete stretchNodes[i];
      delete directionNodes[i];
    }
  }
  else if (min==1)
  {
    for (i = 0; i<Ns; i++)
    {
      delete directionNodes[i];
    }
  }

  return (0);
  
}


