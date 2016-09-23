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

using namespace voom;
using namespace std;

int main(int argc, char* argv[])
{
  time_t start,end;
  time (&start);
  ifstream inp;

  if(argc < 2){
    cout << "Input file missing." << endl;
    return(0);
  }
 
  // File names
  string parameterFileName = argv[1];
  string modelName;
  string OutlineName;
  
  // Minimization constants
  int min = 0;            // Minimize with respect to eta and theta flag (3 -> both, 2 -> only theta, 1-> only eta, 0-> none)
  int Ns_eta = 1;         // Number of stretch magnitude dof
  int Ns_theta = 1;       // Number of stretch direction dof
  double in_eta = 1.0;    // eta initial guess
  double in_theta = 0.0;  // theta initial guess
  double in_phi = 0.0;    // phi initial value
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  // Stretching
  double mu = 1.0;
  double kS = 1.0;
  // Additional energy input parameters
  int WcType = 1;
  double WcConst0 = 1.0;
  double WcConst1 = 1.0;
  double WcConst2 = 1.0;
  double WcConst3 = 1.0;
  double WcConst[4];
  int WcStepsEta = 1;
  int WcStepsTheta = 1;

  // Output
  unsigned int NPout = 0; // Nodes needed to build the outline
  double HexonShift = 0.0;// To separate direction arrows
  string OutFile;
  int refinement;

  // Parameters for MonteCarlo solver
  unsigned int maxIterMTS = 100;
  double finalRatioMTS = 1.0e-5;
  double T1 = 300.0;
  double T2 = 1.0;
  int VarEta = -2;
  int VarTheta = 0;



  bool ferroMagnetic = true;
  vector<vector<unsigned int > > T9ThreeFold(20, vector<unsigned int >(3, 0));
  inp.open("HexamerMapping.dat", ios::in);
  if (!inp) {
    cout << "Cannot open input file HexamersMapping.dat "<< endl;
    return(0);
  }
  for (unsigned int T=0; T<20; T++)
  {
    inp >> T9ThreeFold[T][0];
    T9ThreeFold[T][0] -= 1;
    inp >> T9ThreeFold[T][1];
    T9ThreeFold[T][1] -= 1;
    inp >> T9ThreeFold[T][2];
    T9ThreeFold[T][2] -= 1;
  }
  inp.close();

 



  // Reading input from file passed as argument
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  string dump;
 
  inp >> dump >> modelName;
  inp >> dump >> OutlineName;
  inp >> dump >> min;
  inp >> dump >> Ns_eta;
  inp >> dump >> Ns_theta;
  inp >> dump >> in_eta;
  inp >> dump >> in_theta;
  inp >> dump >> in_phi;
  inp >> dump >> KC;
  inp >> dump >> KG;
  inp >> dump >> C0;
  inp >> dump >> mu;
  inp >> dump >> kS;
  inp >> dump >> WcType;
  inp >> dump >> WcConst0;
  inp >> dump >> WcConst1;
  inp >> dump >> WcConst2;
  inp >> dump >> WcConst3;
  inp >> dump >> NPout;
  inp >> dump >> HexonShift;
  inp >> dump >> OutFile;
  inp >> dump >> refinement;
  inp >> dump >> maxIterMTS;
  inp >> dump >> finalRatioMTS;
  inp >> dump >> T1;
  inp >> dump >> T2;
  inp >> dump >> VarEta;
  inp >> dump >> VarTheta;

  inp.close();

  // Populate Wconst
  WcConst[0] = WcConst0;
  WcConst[1] = WcConst1;
  WcConst[2] = WcConst2;
  WcConst[3] = WcConst3;

  // List input parameters
  cout << " modelName    : " << modelName     << endl
       << " OutlineName  : " << OutlineName   << endl
       << " min          : " << min           << endl
       << " Ns_eta       : " << Ns_eta        << endl
       << " Ns_theta     : " << Ns_theta      << endl
       << " in_eta       : " << in_eta        << endl
       << " in_theta     : " << in_theta      << endl
       << " in_phi       : " << in_phi        << endl
       << " KC           : " << KC            << endl
       << " KG           : " << KG            << endl
       << " C0           : " << C0            << endl
       << " mu           : " << mu            << endl
       << " kS           : " << kS            << endl
       << " WcType       : " << WcType        << endl
       << " alpha        : " << WcConst[0]    << endl
       << " beta         : " << WcConst[1]    << endl
       << " gamma        : " << WcConst[2]    << endl 
       << " delta        : " << WcConst[3]    << endl 
       << " NPout        : " << NPout         << endl
       << " HexonShift   : " << HexonShift    << endl
       << " OutFile      : " << OutFile       << endl
       << " refinement   : " << refinement    << endl
       << " MaxIter      : " << maxIterMTS    << endl
       << " FinalRatio   : " << finalRatioMTS << endl
       << " T1           : " << T1            << endl
       << " T2           : " << T2            << endl
       << " VarEta       : " << VarEta        << endl
       << " VarTheta     : " << VarTheta      << endl; 
 
 


  // vector<double > AngleShift(3, 0.0), ReferenceAngle(3, 0.0); // Needed for T16 from Prof. Rudnick
  // AngleShift[0] = 0.0; AngleShift[1] = M_PI/3.0; AngleShift[2] = 0.0;
  // AngleShift[0] = 2.0*M_PI/3.0; AngleShift[1] = 2.0*M_PI/3.0; AngleShift[2] = M_PI/3.0;
  // To print stretch directions with respect to minimum symmetric configuration
  // ReferenceAngle[0] = 2.0*M_PI/3.0;  ReferenceAngle[1] = 2.0*M_PI/3.0;  ReferenceAngle[2] = M_PI/3.0; // Needed for T16 from Prof. Rudnick
  // ReferenceAngle[0] = 0.0;  ReferenceAngle[1] = 0.0;  ReferenceAngle[2] = 0.0;



 
  // Read in the mesh, in (ascii) legacy vtk format.
  // Create input stream
  ifstream ifs;
  ifs.open(modelName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << modelName << endl;
    exit(0);
  }
 
  // Create vector of nodes
  unsigned int dof = 0, npts = 0, i = 0, j = 0;
  vector<NodeBase* > nodes;
  vector<DeformationNode<3>* > defNodes;
  double Ravg = 0;
  
  // Input .vtk file containing nodes and connectivities
  string token;
  while( token != "POINTS" ) ifs >> token;
  ifs >> npts; 
  defNodes.reserve(npts);
  ifs >> token;   // skip number type
  // read in points
  for(i = 0; i < npts; i++) {
    DeformationNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);

    for(j = 0; j < 3; j++) idx[j] = dof++;

    DeformationNode<3>* n = new DeformationNode<3>(i,idx,x);
    nodes.push_back( n );
    defNodes.push_back( n );
  }
  Ravg /= nodes.size();

  // Read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  std::vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int,3> ct;
  int ntri = 0, tmp = 0, a = 0;
  ifs >> ntri;
  connectivities.reserve(ntri);
  
  ifs >> tmp;
  for (i = 0; i < ntri; i++)
  {
    ifs >> tmp;
    if(tmp != 3) cout << "Some mistake reading the elements connectivity from file. Check again." << endl;
    for(a = 0; a < 3; a++) ifs >> ct[a];
    connectivities.push_back(ct);
  }
  cout << "Number of nodes    = " << nodes.size()          << endl
       << "Ravg               = " << Ravg                  << endl;
  cout << "Number of elements = " << connectivities.size() << endl;
  
  

  
  /*
  // Only for Prof. Rudnick T16
  // Read in element type
  while( token != "eletype" ) ifs >> token;
  std::vector<int > eletype;
  eletype.reserve(ntri);
  ifs >> token;
  ifs >> token;
  ifs >> token;
  cout << token << endl;
  for (i = 0; i < ntri; i++)
  {
    ifs >> tmp;
    eletype.push_back(tmp);
  }
  */



  // The capsomer part goes here with files from Ichosahedron Matlab script
  // Read capsomers to assign different stretch dof to each one
  // Assume that hexons are written before pentons in the vtk input file.
  vector<unsigned int > CapsomersNum(ntri, 0.0);
  while( token != "capsomer") ifs >> token;
  ifs >> token;
  ifs >> token;
  ifs >> token;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >> CapsomersNum[i];
  }
  cout << "Number of capsomers = " <<  CapsomersNum[ntri-1]+1 << endl;



  // Read in the stretch directions. These are vectors pointing from vertex to vertex in each hexon.
  // This vector should point roughly along one of the edges of each triangular element.
  // Since we will define a "canonical" equilateral triangle reference configuration for each element, 
  // we need to find which edge the vector is closest to pointing along, and then choose that edge to
  // define the "stretch angle" defining the reference pre-stretch.
  
  while( token != "shear_direction") ifs >> token; 
  double AvgEdgeLength = 0.0, stretch = 0.0, TOL = 2.5e-1, stretch_dir_norm = 0.0, factor = 0.0;
  vector<double> stretch_angle;
  vector< std::vector<int> > pent_connectivities, hex_connectivities;
  std::vector<int> cm(3,0.0);
  tvmet::Vector<double,3> stretch_dir;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >> stretch_dir(0) >> stretch_dir(1) >> stretch_dir(2);
    stretch_dir_norm = tvmet::norm2(stretch_dir); 
    
    for(j = 0; j < 3; j++) {
      cm[j] = connectivities[i](j);
    }
    
    // Edge vectors in current config.
    tvmet::Vector<double,3> e31(defNodes[cm[0]]->point()-defNodes[cm[2]]->point()), 
                            e32(defNodes[cm[1]]->point()-defNodes[cm[2]]->point()),
                            e12(defNodes[cm[1]]->point()-defNodes[cm[0]]->point()),
                            eCent(defNodes[cm[2]]->point());
    // Compute averate edge length for each triangle
    AvgEdgeLength += (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;

    // pentamers don't have stretch directions
    if(stretch_dir_norm == 0) 
    {
      pent_connectivities.push_back(cm);
    }
    // hexamers do
    else 
    {
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
      // else {std::cout<<"Problem. The stretch direction is not parallel to any triangle edge"<<std::endl; exit(1);}

      // Check that all the normals are defined in the same direction (inner or outward with respect to the capsid surface...
      // otherwise the dof theta does not make sense)
      tvmet::Vector<double,3> n = tvmet::cross(e31,e32);
      if (i==0) {
	factor = (tvmet::dot(n,eCent) > 0) - (tvmet::dot(n,eCent) < 0);
	cout << "tvmet::dot(n,eCent) = " << factor << endl;
      }
	
      if (factor*tvmet::dot(n,eCent) < 0.0 )
      {
	cout << "Theta dof are NOT consistently defined ... ABORT :( " << endl;
	exit(0);
      }	 
      
      // Convert from degrees to radians.
      // stretch_angle.push_back(stretch*M_PI/180.0); // With T16 from Prof. Rudnick
      stretch_angle.push_back(0.0); // With files from Ichosahedron Matlab script
    }
  }
  AvgEdgeLength /= ntri;
  cout << "Number of penton elements = " <<pent_connectivities.size() << endl
       << "Number of hexon elements = "  <<hex_connectivities.size() << endl;

  /* Only with files from Prof. Rudnick 
  // Read capsomers to assign different stretch dof to each one
  // Assume that hexons are written before pentons in the vtk input file.
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
  */

  // Close mesh file in vtk format
  ifs.close();
  




  // Initialize stretch node (ScalarFieldNode type)
  vector<ScalarFieldNode<3>* > stretchNodes, directionNodes;
  stretchNodes.reserve(Ns_eta);   directionNodes.reserve(Ns_theta);
  
  NodeBase::DofIndexMap idE(1);
  ScalarFieldNode<3>::PositionVector pE(0.0);
  unsigned int id = nodes.size();
  // min = 3 -> both, 2 -> only theta, 1-> only eta, 0-> none
  ScalarFieldNode<3>* phiNode = NULL;
  if (min == 0 || min == 1 || min ==3)
  {
    for (i = 0; i < Ns_eta; i++)
    {
      idE[0] = dof++;
      ScalarFieldNode<3>* stretchN = new ScalarFieldNode<3>(id, idE, pE, in_eta);
      id++;
      stretchNodes.push_back(stretchN);
    }
    // ! To keep some of the dof constant we can decide to not insert them in Model BUT
    // THIS ONLY WORKS IF the constant dof are numbered last
    for (i = 0; i < Ns_theta; i++)
    {
      idE[0] = dof++;
      ScalarFieldNode<3>* dirN = new ScalarFieldNode<3>(id, idE, pE, in_theta);
      id++;
      directionNodes.push_back(dirN);
    }
    if (WcType == 5) {
      idE[0] = dof++;
      phiNode = new ScalarFieldNode<3>(id, idE, pE, in_phi);
      id++;
    }
  }
  else
  {
    // ! To keep some of the dof constant we can decide to not insert them in Model BUT
    // THIS ONLY WORKS IF the constant dof are numbered last
    for (i = 0; i < Ns_theta; i++)
    {
      idE[0] = dof++;
      ScalarFieldNode<3>* dirN = new ScalarFieldNode<3>(id, idE, pE, in_theta);
      id++;
      directionNodes.push_back(dirN);
    }
    if (WcType == 5) {
       idE[0] = dof++;
       phiNode = new ScalarFieldNode<3>(id, idE, pE, in_phi);
       id++;
    }
    for (i = 0; i < Ns_eta; i++)
    {
      idE[0] = dof++;
      ScalarFieldNode<3>* stretchN = new ScalarFieldNode<3>(id, idE, pE, in_eta);
      id++;
      stretchNodes.push_back(stretchN);
    }
  }  


 

  // Create bending body; SCElastic only has curvature terms, no in-plane strain terms
  int quadOrder = 1;
  // Bending material object, to be copied when body generates new elements
  typedef SCElastic BendingMaterial;
  BendingMaterial bending(KC,KG,C0);

  // Bending body will generate a mesh of Loop Subdivision shell elements from the connectivity and nodes
  typedef LoopShellBody<BendingMaterial> LSB;
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);
  // Use the body's output routines to generate vtk files with results
  bending_body->setOutput(paraview);
  bending_body->compute(true, false, false);
  cout << "En_Bend = " << bending_body->totalStrainEnergy() << endl;
  cout << "En_Bend = " << bending_body->energy() << endl;
  // bending_body->checkConsistency();

  

  // Create two stretching bodies, one for pents and one for hex's, which don't have any bending energy
  // (bending moduli set to zero)
  // Use EvansElastic Material, which separates area strain from pure stretch strain (at constant area).  


  // Pentamers
  typedef EvansElastic PentStretchingMaterial;
  PentStretchingMaterial PentStretch(0.0, 0.0, 0.0, mu, kS);

  // Use generic C0MembraneBody for stretching of pents, with the EvansElastic material.
  typedef C0MembraneBody<TriangleQuadrature, PentStretchingMaterial, ShapeTri3> pent_membrane;
  pent_membrane * pentamer_body = new pent_membrane(PentStretch, pent_connectivities, nodes, quadOrder);
  pentamer_body->setOutput(paraview);
  pentamer_body->compute(true, false, false);
  cout << "En_Pent = " << pentamer_body->totalStrainEnergy() << endl;
  cout << "En_Pent = " << pentamer_body->energy() << endl;
  // pentamer_body->checkConsistency();
  
 
  
  // Hexamers
  // Use specialezed StretchHexonBody (modified version of C0MembraneBody), which has material hard-coded as
  // ModifiedEvansElasticMin, which takes care of the multiplicative F=AG "Eigen-strain" decomposition.



  // Shift in initial hexamer configuration
  inp.open("T9ShiftAnglesInfCur.dat", ios::in);
  if (!inp) {
    cout << "Cannot open input file AngleShiftFile" << endl;
    return(0);
  }
  
  // For T9
  vector<double > VAngleShift(stretchNodes.size(), 0.0);
  for (i = 0; i < VAngleShift.size(); i++)
  {
    inp >> VAngleShift[i];
    VAngleShift[i] *= (M_PI/3.0);
    // cout << VAngleShift[i] << endl;
  }
  
  
  unsigned int CurrentCap = 0;
  unsigned int ind = 0;
  for (unsigned int elT = 0; elT < stretch_angle.size(); elT++)
  {
    if (CapsomersNum[elT] != CurrentCap) {
      ind++;
      CurrentCap = CapsomersNum[elT];
      // cout << CurrentCap << endl;
    }
    stretch_angle[elT] = VAngleShift[ind]; 
    // cout <<  stretch_angle[elT] << endl;
  }



  // Correct the stretch modulus if using double-well potential.
  // Not necessary if using single-well.
  // mu = mu/(2*(eta*eta+c2)+2*c1);
  TriangleQuadrature  HexQuad(1);
  const ShapeTri3::CoordinateArray dummy(0.0);
  ShapeTri3 HexShape(dummy);
  StretchHexonBody * hex_body = new StretchHexonBody(hex_connectivities, defNodes, stretchNodes, directionNodes, stretch_angle, mu, kS, WcType, WcConst, HexQuad, &HexShape, CapsomersNum, phiNode);
  hex_body->setOutput(paraview);
  hex_body->compute(true, false, false);
  cout << "En_Hex = " << hex_body->totalStrainEnergy() << endl;
  cout << "En_Hex = " << hex_body->energy() << endl;
  cout << "En_Hex = " << hex_body->totalConformationalEnergy() << endl;
  // hex_body->checkConsistency();
  

  
  // Set reference configuration
  // Define the edge length of the canonical 
  // const double EDGELEN = 2.0; - > AvgEdgeLength is used instead
  
  bending_body->SetRefConfiguration(AvgEdgeLength);
  pentamer_body->SetRefConfiguration(AvgEdgeLength);
  hex_body->SetRefConfiguration(AvgEdgeLength);
   
  //   pentamer_body->SetRefConfiguration(2.0);
  //   hex_body->SetRefConfiguration(2.0);




  switch (min) {
  case 1: // stretch magnitude is a dof
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    break;
  case 2: // stretch direction is a dof
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
    if (WcType == 5) {
      nodes.push_back(phiNode);
    }
    break;
  case 3: // stretch direction is a dof
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
    if (WcType == 5) {
      nodes.push_back(phiNode);
    }
    break;
  }
 
  // Hard coded for T12
  vector<double > VRefAngle(110, 0.0), ZeroAngles(110, 0.0);
  double angleRef = M_PI/3.0;
  VRefAngle[0] = 0.0;  VRefAngle[1] = 0.0;  VRefAngle[2] = angleRef;  VRefAngle[3] = 0.0;  VRefAngle[4] = angleRef;
  VRefAngle[5] = angleRef;  VRefAngle[6] = 2.0*angleRef;  VRefAngle[7] = 0.0;  VRefAngle[8] = angleRef;  VRefAngle[9] = 0.0;
  VRefAngle[10] = angleRef;  VRefAngle[11] = 2.0*angleRef;  VRefAngle[12] = 2.0*angleRef;  VRefAngle[13] =-angleRef;  VRefAngle[14] = 0.0;
  VRefAngle[15] = angleRef;  VRefAngle[16] = 0.0;  VRefAngle[17] = 2.0*angleRef;  VRefAngle[18] = 2.0*angleRef;  VRefAngle[19] = -angleRef;
  
  VRefAngle[20] = 0.0;  VRefAngle[21] = angleRef;  VRefAngle[22] = 0.0;  VRefAngle[23] = angleRef;  VRefAngle[24] = 2.0*angleRef;
  VRefAngle[25] =-angleRef;  VRefAngle[26] = 0.0;  VRefAngle[27] = angleRef;  VRefAngle[28] = 0.0;  VRefAngle[29] = 2.0*angleRef;
  VRefAngle[30] = 0.0;  VRefAngle[31] = angleRef;  VRefAngle[32] = 0.0;  VRefAngle[33] = angleRef;  VRefAngle[34] = 2.0*angleRef;
  VRefAngle[35] =2.0*angleRef;  VRefAngle[36] = 0.0;  VRefAngle[37] = angleRef;  VRefAngle[38] = 0.0;  VRefAngle[39] = 0.0;
  
  VRefAngle[40] = angleRef;  VRefAngle[41] = 2.0*angleRef;  VRefAngle[42] = 0.0;  VRefAngle[43] = angleRef;  VRefAngle[44] = 0.0;
  VRefAngle[45] = angleRef;  VRefAngle[46] = 2.0*angleRef;  VRefAngle[47] = 0.0;  VRefAngle[48] = angleRef;  VRefAngle[49] = 0.0;
  VRefAngle[50] = 0.0;  VRefAngle[51] = angleRef;  VRefAngle[52] = 2.0*angleRef;  VRefAngle[53] = 0.0;  VRefAngle[54] = angleRef;
  VRefAngle[55] = 0.0;  VRefAngle[56] = 0.0;  VRefAngle[57] = 2.0*angleRef;  VRefAngle[58] = -angleRef;  VRefAngle[59] = 0.0;

  VRefAngle[60] = angleRef;  VRefAngle[61] = 0.0;  VRefAngle[62] = 2.0*angleRef;  VRefAngle[63] = 2.0*angleRef;  VRefAngle[64] = 0.0;
  VRefAngle[65] = angleRef;  VRefAngle[66] = 0.0;  VRefAngle[67] = angleRef;  VRefAngle[68] = angleRef;  VRefAngle[69] = 2.0*angleRef;
  VRefAngle[70] = 0.0;  VRefAngle[71] = angleRef;  VRefAngle[72] = 0.0;  VRefAngle[73] = angleRef;  VRefAngle[74] = 2.0*angleRef;
  VRefAngle[75] = 0.0;  VRefAngle[76] = angleRef;  VRefAngle[77] = 0.0;  VRefAngle[78] = angleRef;  VRefAngle[79] = 2.0*angleRef;

  VRefAngle[80] =-angleRef;  VRefAngle[81] = 0.0;  VRefAngle[82] = angleRef;  VRefAngle[83] = 0.0;  VRefAngle[84] = angleRef;
  VRefAngle[85] = 2.0*angleRef;  VRefAngle[86] = 0.0;  VRefAngle[87] = angleRef;  VRefAngle[88] = 0.0;  VRefAngle[89] = 2.0*angleRef;
  VRefAngle[90] = 2.0*angleRef;  VRefAngle[91] = 0.0;  VRefAngle[92] = 0.0;  VRefAngle[93] = angleRef;  VRefAngle[94] = 0.0;
  VRefAngle[95] = angleRef;  VRefAngle[96] = 2.0*angleRef;  VRefAngle[97] = 0.0;  VRefAngle[98] = angleRef;  VRefAngle[99] = 0.0;

  VRefAngle[100] = 2.0*angleRef;  VRefAngle[101] = 0.0;  VRefAngle[102] = 0.0;  VRefAngle[103] = angleRef;  VRefAngle[104] = 0.0;
  VRefAngle[105] = 2.0*angleRef;  VRefAngle[106] = 0.0;  VRefAngle[107] = angleRef;  VRefAngle[108] = 0.0;  VRefAngle[109] = 2.0*angleRef;


  /* Needed for T16 from Prof. Rudnick
  // compute new angle offset - Number of examers is hard coded !!
  // For T16
      vector<double > VectorReferenceAngle(150, 0.0), VectorAngleShift(150, 0.0);
      for (unsigned int exT = 0; exT < 150; exT++)
	{
	  if (exT < 60) {
	    VectorReferenceAngle[exT] = ReferenceAngle[0]; 
	    VectorAngleShift[exT] = AngleShift[0]; 
	  }
	  else if (exT< 90) {
	    VectorReferenceAngle[exT] = ReferenceAngle[1];
	    VectorAngleShift[exT] = AngleShift[1];  
	  }
	  else {
	    VectorReferenceAngle[exT] = ReferenceAngle[2]; 
	    VectorAngleShift[exT] = AngleShift[2]; 
	  }
	}
      vector<double > AngleOffset(stretch_angle.size(), 0.0);
      for (unsigned int elT = 0; elT < stretch_angle.size(); elT++)
	{
	  if (eletype[elT] == 61)
	    {
	      AngleOffset[elT] = stretch_angle[elT] + AngleShift[0];
	    }
	  else if (eletype[elT] == 62)
	    {
	      AngleOffset[elT] = stretch_angle[elT] + AngleShift[1];
	    }
	  else
	    {
	      AngleOffset[elT] = stretch_angle[elT] + AngleShift[2];
	    }
	}
      
  // assign new angle offset
  hex_body->setAngleOffset(AngleOffset);
  */




  // Create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(pentamer_body);
  bdc.push_back(hex_body);
 
  Model model(bdc, nodes);
  // Consistency check
  // model.checkConsistency(true,false);
  // model.checkRank(model.dof()-6,true);


  
  
  // Set solver and solve
  int iprint = 0;
  double factr = 1.0e1;
  double pgtol = 1.0e-6;
  int m = 5;
  int maxIter = 20000;
  cout << "Input iprint: " << iprint << endl
       << "Input factr:  " << factr  << endl
       << "Input pgtol:  " << pgtol  << endl
       << "Input m:      " << m      << endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);

  PrintingStretches PrintTvirus(modelName, OutlineName, NPout, 
				OutFile+"iter", OutFile+"outline", OutFile+"directions", OutFile+"elemdir", 1,
			   	bdc, stretchNodes, directionNodes, CapsomersNum, 
				VAngleShift, ZeroAngles, HexonShift, WcConst[1], AvgEdgeLength, refinement, false, T9ThreeFold, 0);

   vector<double > SymmIndex = PrintTvirus.printMaster(-1);
   cout << "SymmIndex [nodes, eta, theta] = " << SymmIndex[0] << " " << SymmIndex[1] << " " << SymmIndex[2] << endl;
   PrintTvirus.setPrintingConfig(true);
   

   // Compute FvK number before solving
   double FvK = 0.0, Y = 0.0;
   Y = 4.0*kS*mu/(kS+mu);
   FvK = Y*Ravg*Ravg/KC;
   cout << " FvK number before solving = " << FvK           << endl;

 


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
  MontecarloTwoStages MTS(MontecarloDof, VarType, &solver, &PrintTvirus, maxIterMTS, false, ferroMagnetic); 
  // MTS.SetTempSchedule(MTS.EXPONENTIAL, T1, T2, finalRatioMTS);
  MTS.SetTempSchedule(MTS.CONSTANT, T1, T2, finalRatioMTS);
  MTS.solve(&model); 
  

  

  // Printing results
  // Print the postprocessing strains and stress
  // hex_body->printParaviewEigVec("Hexbody");
  // pentamer_body->printParaviewEigVec("Pentbody");
  // hex_body->printParaview2("Hexbody2");
  // pentamer_body->printParaview2("Pentbody2");

 


  // Calculate capsid center, average radius and asphericity
  tvmet::Vector<double,3> center(0.0);
  Ravg = 0.0;
  double deltaR = 0.0, Rtemp = 0.0, asphericity = 0.0;
  for(i = 0; i < defNodes.size(); i++)
  {
    center += defNodes[i]->point();
    Ravg += tvmet::norm2(defNodes[i]->point());
  }
  center /= defNodes.size();
  Ravg /= defNodes.size();
  cout << "Center      = " << center  << endl;
  

  for(i = 0; i < defNodes.size(); i++)
  {
    Rtemp = tvmet::norm2(defNodes[i]->point());
    deltaR += pow(Rtemp-Ravg,2);
  }
  asphericity = sqrt( deltaR/defNodes.size() )/Ravg;


  // Compute max difference in eta and theta
  double EtaMaxDiff = 0.0, EtaTempDiff = 0.0;
  if (Ns_eta > 1)
  {
    EtaMaxDiff = abs(stretchNodes[1]->point()-stretchNodes[0]->point());
    for (i = 1; i < Ns_eta-1; i++)
    {
      EtaTempDiff = abs(stretchNodes[i+1]->point()-stretchNodes[i]->point());
      if ( EtaTempDiff > EtaMaxDiff ) EtaMaxDiff = EtaTempDiff;
    }
  }




  //---------------//
  // Print results //
  // Print the total energies of three bodies
  
  double Wben = 0.0, Wpen = 0.0, Whex = 0.0, Wtot = 0.0, Wconf = 0.0;
  FvK = Y*Ravg*Ravg/KC;
  Wben  =  bending_body->totalStrainEnergy();
  Wpen  =  pentamer_body->totalStrainEnergy();
  Whex  =  hex_body->totalStrainEnergy();
  Wconf =  hex_body->totalConformationalEnergy();
  Wtot = Wben + Wpen + Whex;
  cout << "----------------------------------------" << endl;

  cout << " Stretch magnitude     = " << endl;
  for (i = 0; i < Ns_eta; i++)
  {
    cout << stretchNodes[i]->point() << "  ";
    if ( (i+1)%5 == 0 || i == Ns_eta-1 ) cout << endl;
  }

  cout << " Direction angle     = " << endl;
  for (i = 0; i < Ns_theta; i++)
  {
    cout << directionNodes[i]->point() << "  ";
    if ( (i+1)%5 == 0 || i == Ns_theta-1 ) cout << endl;
  }
  
  cout << " Max diff. in eta    = " << EtaMaxDiff    << endl;
  cout << " FvK number          = " << FvK           << endl;
  cout << " Bending energy      = " << Wben          << endl;
  cout << " Pent Stretch energy = " << Wpen          << endl;
  cout << " Hex  Stretch energy = " << Whex          << endl;
  cout << " Hex conform. energy = " << Wconf         << endl;
  cout << " Total energy        = " << Wtot          << endl;
  cout << " Total def. energy   = " << Wtot - Wconf  << endl;
  cout << " Asphericity         = " << asphericity   << endl;
  cout << " Ravg                = " << Ravg          << endl;
  cout << " AverageEdgeLength   = " << AvgEdgeLength << endl;
  cout << "----------------------------------------" << endl;
  
  time (&end);
  double dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;
 

  
  // Clean up
  delete bending_body;
  delete pentamer_body;
  delete hex_body;

  for (i = 0; i<nodes.size(); i++)
  {
    delete nodes[i];
  }

  switch (min) {
  case 0: 
    for (i = 0; i<Ns_eta; i++)
    {
      delete stretchNodes[i];
    }
    for (i = 0; i<Ns_theta; i++)
    {
      delete directionNodes[i];
    }
    break;
  case 1:
    for (i = 0; i<Ns_theta; i++)
    {
      delete directionNodes[i];
    }    
    break;
  case 2:
    for (i = 0; i<Ns_eta; i++)
    {
      delete stretchNodes[i];
    }
    break;
  }

  return (0);
  
}







    
    
   
