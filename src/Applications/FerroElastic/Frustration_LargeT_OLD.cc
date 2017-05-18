#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <tvmet/Vector.h>
#include <cmath>
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
#include "VoomMath.h"
#include "Utils/PrintingStretches.h"

using namespace voom;
using namespace std;

/*! This program creates a shell model with LoopShellBody for bending
  and 2 C0MembraneBodies for stretching, one of pentamers and one of
  hexamers.
 */
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
  int min = 0;            // Minimize with respect to eta and theta flag (3 -> both, 2 -> only theta, 1-> only eta, 0-> none)
  int Ns_eta = 1;         // Number of stretch magnitude dof
  int Ns_theta = 1;       // Number of stretch direction dof
  double in_eta = 0.0;    // eta initial guess
  double in_theta = 0.0;  // theta initial guess
  double in_phi = 0.0;    // phi initial value
 
  // Bending
  double KC = 1.0;
  double KG =-1.0;
  double C0 = 0.0;
  // Stretching
  double mu = 1.0;
  double kS = 1.0;
   
  // Output
  unsigned int NPout = 0; // Nodes needed to build the outline
  string OutFile;
  int refinement = 1;

  // Input file for refinement angle
  string AngleShiftFile;

  // Analysis Type
  int AnType = 0;      // 0 -> Shield analyses, 1 -> Energy analyses (three fold sites hexamers in T9), 2 -> Energy analyses (all hexamers)
  int NumEnergy = 1;  // number of energies to be computed for different three fold sites configurations

  // Reading input from file passed as argument
  ifstream inp;
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
  inp >> dump >> NPout;
  inp >> dump >> OutFile;
  inp >> dump >> refinement;
  inp >> dump >> AngleShiftFile;
  inp >> dump >> AnType;
  inp >> dump >> NumEnergy;
  
  inp.close();

  // Populate WcType, WcConst, HexonShift (not used in these analyses)
  int WcType = 0;
  double WcConst[4];
  WcConst[0] = 0.0;
  WcConst[1] = 0.0;
  WcConst[2] = 0.0;
  WcConst[3] = 0.0;
  double HexonShift = 0.1;

  // List input parameters
  cout << " modelName      : " << modelName      << endl
       << " OutlineName    : " << OutlineName    << endl
       << " min            : " << min            << endl
       << " Ns_eta         : " << Ns_eta         << endl
       << " Ns_theta       : " << Ns_theta       << endl
       << " in_eta         : " << in_eta         << endl
       << " in_theta       : " << in_theta       << endl
       << " in_phi         : " << in_phi         << endl
       << " KC             : " << KC             << endl
       << " KG             : " << KG             << endl
       << " C0             : " << C0             << endl
       << " mu             : " << mu             << endl
       << " kS             : " << kS             << endl
       << " NPout          : " << NPout          << endl
       << " OutFile        : " << OutFile        << endl
       << " Refinement     : " << refinement     << endl
       << " AngleShiftFile : " << AngleShiftFile << endl
       << " AnType         : " << AnType         << endl
       << " NumEnergy      : " << NumEnergy      << endl;
      
  



  // Analyses list
  vector<unsigned int> NearNeighList(4, 0);
  NearNeighList[0] = 3;    // Three fold site (TFS)
  NearNeighList[1] = 8;    // Surrounding three fold site
  NearNeighList[2] = 24;   // Surrounding three fold site
  NearNeighList[3] = 75;   // Surrounding three fold site
  
  vector<vector<unsigned int > > TestShift(81, vector<unsigned int >(4,0));
  unsigned int i = 0, j = 0, k = 0, l = 0, ind = 0;
  for (i = 0; i<3; i++) {
    for (j = 0; j<3; j++) {
      for (k = 0; k<3; k++) {
	for (l = 0; l<3; l++) {
	  TestShift[ind][0] = l;
	  TestShift[ind][1] = k;
	  TestShift[ind][2] = j;
	  TestShift[ind][3] = i;
	  ind++;
	}
      }
    }
  }




  // Read in the mesh, in (ascii) legacy vtk format.
  // Create input stream
  ifstream ifs;
  ifs.open(modelName.c_str(), ios::in);
  if (!ifs) {
    cout << "Cannot open input file: " << modelName << endl;
    exit(0);
  }
 
  // Create vector of nodes
  vector<NodeBase* > nodes;
  vector<DeformationNode<3>* > defNodes;
  double Ravg = 0;
  unsigned int dof = 0, npts = 0;
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
  vector<double > ZeroAngles;
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
      /*else {
	std::cout<<"Problem. The stretch direction is not parallel to any triangle edge"<<std::endl; exit(1);
	}*/

      // Check that all the normals are defined in the same direction (inner or outward with respect to the capsid surface...
      // otherwise the theta dof does not make sense)
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
      // stretch_angle.push_back((stretch)*M_PI/180.); 
      // cout << stretch << endl;
      ZeroAngles.push_back(0.0);
      stretch_angle.push_back(0.0); 
    }
  }
  AvgEdgeLength /= ntri;
  cout << "Number of penton elements = " <<pent_connectivities.size() << endl
       << "Number of hexon elements = "  <<hex_connectivities.size() << endl;

  
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
  // blitz::Array<double,1> LowerBound, UpperBound;
  // blitz::Array<int,1> BoundType;
  int NumDOF = 0;
  // Determine number of DOF we will solve for
  switch (min)
  {
    case 0:
      NumDOF = nodes.size()*3;
      break;
    case 1:
      NumDOF = nodes.size()*3 + Ns_eta;
      break;
    case 2:
      NumDOF = nodes.size()*3 + Ns_theta;
      if (WcType == 5 || WcType == 6) {
	NumDOF += 1; // including phi dof
      }
      break;
    case 3:
      NumDOF = nodes.size()*3 + Ns_eta + Ns_theta;
      if (WcType == 5 || WcType == 6) {
	NumDOF += 1; // including phi dof
      }
      break;
    default:
      cout << "Minimization case not implemented!" << endl;
      return(0);
  }



  if (min == 0 || min == 1 || min ==3)
  {
    // Initialize bound of dof
    // LowerBound.resize(NumDOF);
    // UpperBound.resize(NumDOF);
    // BoundType.resize(NumDOF);
    // LowerBound = 0.0;
    // UpperBound = 0.0;
    // BoundType  = 0;  // Default is no bound

    for (i = 0; i < Ns_eta; i++)
    {
      idE[0] = dof++;
      ScalarFieldNode<3>* stretchN = new ScalarFieldNode<3>(id, idE, pE, in_eta);
      // if (min != 0) {
	// BoundType(idE[0]) = 1;  // Only lower bound
	// LowerBound(idE[0]) = 1.0 + WcConst[1];
      // }
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
    if (WcType == 5 || WcType == 6) {
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
    if (WcType == 5 || WcType == 6) {
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
  // ModifiedEvansElastic, which takes care of the multiplicative F=AG "Eigen-strain" decomposition.

  // Correct the stretch modulus if using double-well potential.
  // Not necessary if using single-well.
  // mu = mu/(2*(eta*eta+c2)+2*c1);
  TriangleQuadrature  HexQuad(1);
  const ShapeTri3::CoordinateArray dummy(0.0);
  ShapeTri3 HexShape(dummy);








  // Shift in initial hexamer configuration
  inp.open(AngleShiftFile.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << AngleShiftFile << endl;
    return(0);
  }
  
  // For T9
  vector<double > VAngleShift(directionNodes.size(), 0.0);
  for (i = 0; i < VAngleShift.size(); i++)
  {
    inp >> VAngleShift[i];
    VAngleShift[i] *= (M_PI/3.0);
    // cout << VAngleShift[i] << endl;
  }
  
  
  unsigned int CurrentCap = 0;
  ind = 0;
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
      
 

  StretchHexonBody * hex_body = new StretchHexonBody(hex_connectivities, defNodes, stretchNodes, directionNodes, stretch_angle, mu, kS, WcType, WcConst, HexQuad, &HexShape, CapsomersNum, phiNode);
  hex_body->setOutput(paraview);
  hex_body->compute(true, false, false);
  cout << "En_Hex    = " << hex_body->totalStrainEnergy() << endl;
  cout << "En_Hex    = " << hex_body->energy() << endl;
  cout << "Wconf_Hex = " << hex_body->totalConformationalEnergy() << endl;
  // hex_body->checkConsistency();
  

  
  // Set reference configuration
  // Define the edge length of the canonical 
  // const double EDGELEN = 2.0; - > AvgEdgeLength is used instead
  
  bending_body->SetRefConfiguration(AvgEdgeLength);
  pentamer_body->SetRefConfiguration(AvgEdgeLength);
  hex_body->SetRefConfiguration(AvgEdgeLength);
 
  cout << "Reference configuration has been imposed " << endl;
   
  //   pentamer_body->SetRefConfiguration(2.0);
  //   hex_body->SetRefConfiguration(2.0);




  /*
NearNeighList[3] = 75;   // Surrounding three fold site
  
  vector<vector<unsigned int > > TestShift(81, vector<unsigned int >(4,0));
  unsigned int i = 0, j = 0, k = 0, l = 0, ind = 0;
  for (j = 0; j<3, j++) {
    for (k = 0; k<3, k++) {
      for (l = 0; l<3, l++) {
	for (i = 0; i<3, i++) {
	  TestShift[ind]
  // Shift in initial hexamer configuration
  inp.open(RefAngleFile.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << RefAngleFile << endl;
    return(0);
  }
  
  // For T9
  vector<double > VRefAngle(stretchNodes.size(), 0.0);
  for (i = 0; i < VRefAngle.size(); i++)
  {
    inp >> VRefAngle[i];
    VRefAngle[i] *= (M_PI/3.0);
    cout << VRefAngle[i] << endl;
  }
  
  vector<double > AngleOffset(stretch_angle.size(), 0.0);
  unsigned int CurrentCap = 0, ind = 0;
  for (unsigned int elT = 0; elT < stretch_angle.size(); elT++)
  {
    if (CapsomersNum[elT] != CurrentCap) {
      ind++;
      CurrentCap = CapsomersNum[elT];
      cout << CurrentCap << endl;
    }
    AngleOffset[elT] = VRefAngle[ind]; 
  }
      
  // assign new angle offset
  hex_body->setAngleOffset(AngleOffset);
  */



    
  switch (min) {
  case 1: // stretch magnitude is a dof
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    break;
  case 2: // stretch direction is a dof
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
    if (WcType == 5 || WcType == 6) {
      nodes.push_back(phiNode);
    }
    break;
  case 3: // stretch direction is a dof
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
    if (WcType == 5 || WcType == 6) {
      nodes.push_back(phiNode);
    }
    break;
  }


  
  

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
  double pgtol = 1.0e-8; // Previously 1.0e-6
  int m = 10;            // Previously 5
  int maxIter = 30000;
  std::ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m;
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
  cout << "Bound  DOF = " << NumDOF << endl;
  cout << "Solver DOF = " << model.dof() << endl;


  // Fill in angle vectors needed for output
  vector<double > SymmIndex(6, 0.0);
  vector<vector<unsigned int > > NOTused;
  PrintingStretches PrintTvirus(modelName, OutlineName, NPout, 
				OutFile+"iter", OutFile+"outline", OutFile+"directions", OutFile+"elemdir", 1,
				bdc, stretchNodes, directionNodes, CapsomersNum, 
				VAngleShift, ZeroAngles, HexonShift, WcConst[1], AvgEdgeLength, refinement, false, NOTused, 0);
  SymmIndex = PrintTvirus.printMaster(0);
  cout << "Initial symmetry indices " << endl;
       cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
       cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
       cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
       cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
       cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
       cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;

  PrintTvirus.setPrintingConfig(true);


   // Compute FvK number before solving
   double FvK = 0.0, Y = 0.0;
   Y = 4.0*kS*mu/(kS+mu);
   FvK = Y*Ravg*Ravg/KC;
   cout << " FvK number before solving = " << FvK           << endl;




   srand(time(NULL));
   // Solve to minimize the energy
   if (AnType == 0)  // Frustration case
   {
     
     // Other Three fold site hexamers
     set<unsigned int > OtherTFS;
     OtherTFS.insert(15); OtherTFS.insert(20); OtherTFS.insert(26); OtherTFS.insert(31);
     OtherTFS.insert(35); OtherTFS.insert(39); OtherTFS.insert(43); OtherTFS.insert(49);
     OtherTFS.insert(52); OtherTFS.insert(56); OtherTFS.insert(59); OtherTFS.insert(65);
     OtherTFS.insert(68); OtherTFS.insert(72); OtherTFS.insert(78); OtherTFS.insert(79);

     /*
     // Select 16 random numbers (equivalent to the number of other 3 fold sites)
     map<uint, double> randBackGround;
     
     cout << "Randnumber = ";
     for (set<uint>::iterator r = OtherTFS.begin(); r != OtherTFS.end(); r++) {
       double randnumber = double(rand()%3);
       randBackGround.insert(make_pair(*r, randnumber));
       VAngleShift[*r] = randnumber*M_PI/3.0;
       cout << *r << " " << randnumber << endl;
     }
     cout << endl;
     
     // Change randomly the background of T9
     for (j = 0; j < stretch_angle.size(); j++) {
       if (OtherTFS.find(CapsomersNum[j]) != OtherTFS.end()) 
       {
	 stretch_angle[j] = randBackGround[CapsomersNum[j]]*M_PI/3.0;
       }
     }
     */

     // Set all other 3 fold site to have zero stretch - An_I
     for (j = 0; j < stretchNodes.size(); j++) {
       if (OtherTFS.find(j) != OtherTFS.end()) 
       {
	 stretchNodes[j]->setPoint(1.0);
       }
     }


     
     // Explore 81 cases for one triplet
     for (i = 0; i < TestShift.size(); i++) {
       for (j = 0; j < stretch_angle.size(); j++) {
	 if (CapsomersNum[j] == NearNeighList[0]) 
	 {
	   stretch_angle[j]             = TestShift[i][0]*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = TestShift[i][0]*M_PI/3.0; 
	 }
	 else if (CapsomersNum[j] == NearNeighList[1]) 
	 {
	   stretch_angle[j]             = TestShift[i][1]*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = TestShift[i][1]*M_PI/3.0; 
	 }
	 else if (CapsomersNum[j] == NearNeighList[2]) 
	 {
	   stretch_angle[j]             = TestShift[i][2]*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = TestShift[i][2]*M_PI/3.0; 
	 }
	 else if (CapsomersNum[j] == NearNeighList[3]) 
	 {
	   stretch_angle[j]             = TestShift[i][3]*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = TestShift[i][3]*M_PI/3.0; 
	 }
       }
       
       // assign new angle offset
       hex_body->setAngleOffset(stretch_angle);

       // Solve 
       solver.solve(&model);
       double Wtot = bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy();
       cout << "STATE = " << i << endl;
       cout << "State = " << TestShift[i][0] << " " << TestShift[i][1] << " " << TestShift[i][2] << " " << TestShift[i][3] << " " << Wtot << endl; 
       cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
       
       SymmIndex = PrintTvirus.printMaster(i);
       cout << "Symmetry Index " << endl;
       cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
       cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
       cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
     }



     // Testing far away shielding in 20 cases
         // Setting neighbor three fold sites in the 112 configuration
         // Setting neighbor three fold sites in the 002 configuration
         for (j = 0; j < stretch_angle.size(); j++) {
	   if (CapsomersNum[j] == NearNeighList[1]) 
	   {
	     stretch_angle[j]             = 0.0; // M_PI/3.0; // 0.0;
	     VAngleShift[CapsomersNum[j]] = 0.0; // M_PI/3.0; // 0.0;
	   }
	   else if (CapsomersNum[j] == NearNeighList[2]) 
	   {
	     stretch_angle[j]             = 0.0; // M_PI/3.0; // 0.0;
	     VAngleShift[CapsomersNum[j]] = 0.0; // M_PI/3.0; // 0.0;
	   }
	   else if (CapsomersNum[j] == NearNeighList[3]) 
	   {
	     stretch_angle[j]             = 2.0*M_PI/3.0;
	     VAngleShift[CapsomersNum[j]] = 2.0*M_PI/3.0; 
	   }
	 }

	 /*	 
     // Look at 20 randomly generated cases 
     srand(time(NULL));
     for (i=0; i<20; i++)
     {
       // Change one Three fold site in the background at the time
       /*
       set<unsigned int >::iterator OtherTFSit = OtherTFS.begin();
       advance(OtherTFSit, rand()%16);
       for (j = 0; j < stretch_angle.size(); j++) {
	 if ( *OtherTFSit == CapsomersNum[j] ) 
	 {
	   cout << "Changed TFS " << *OtherTFSit << endl;
	   double randNumber = double(rand()%3);
	   stretch_angle[j]             = randNumber*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = randNumber*M_PI/3.0; 
	 }
       }
       
       
       // Change All three fold sites in the background at once
       cout << "randnumber = ";
       randBackGround.clear();
       for (set<uint>::iterator r = OtherTFS.begin(); r != OtherTFS.end(); r++) {
	 double randnumber = double(rand()%3);
	 randBackGround.insert(make_pair(*r, randnumber));
	 VAngleShift[*r] = randnumber*M_PI/3.0;
	 cout << *r << " " << randnumber << endl;
       }
       cout << endl;
     
       // Change randomly the background of T9
       for (j = 0; j < stretch_angle.size(); j++) {
	 if (OtherTFS.find(CapsomersNum[j]) != OtherTFS.end()) 
	   {
	     stretch_angle[j] = randBackGround[CapsomersNum[j]]*M_PI/3.0;
	   }
       }
       
       




       for (k = 0; k < 3; k++) {
	 for (j = 0; j < stretch_angle.size(); j++) {
	   if (CapsomersNum[j] ==  NearNeighList[0]) 
	   {
	     stretch_angle[j]             = k*M_PI/3.0;
	     VAngleShift[CapsomersNum[j]] = k*M_PI/3.0; 
	   }
	 }
       
	 // assign new angle offset
	 hex_body->setAngleOffset(stretch_angle);

	 // Solve 
	 solver.solve(&model);
	 double Wtot = bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy();
	 cout << "STATE = " << i*3 + k + 81 << endl;
	 // cout << "State = " << k << " " << TestShift[80][1] << " " << TestShift[80][2] << " " << TestShift[80][3] << " " << Wtot << endl; 
	 cout << "State = " << k << " 1 1 2 " << Wtot << endl;
	 cout << "Wconf = " << hex_body->totalConformationalEnergy() << endl;
       
	 SymmIndex = PrintTvirus.printMaster(i*3 + k + 81);
	 cout << "Symmetry Index " << endl;
	 cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	 cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	 cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
       }
     }*/
 
     
   }
   else if(AnType == 1)
   {
     // Three fold sites hexamers
     set<unsigned int > OtherTFS;
     OtherTFS.insert(3);  OtherTFS.insert(8);  OtherTFS.insert(24); OtherTFS.insert(75);
     OtherTFS.insert(15); OtherTFS.insert(20); OtherTFS.insert(26); OtherTFS.insert(31);
     OtherTFS.insert(35); OtherTFS.insert(39); OtherTFS.insert(43); OtherTFS.insert(49);
     OtherTFS.insert(52); OtherTFS.insert(56); OtherTFS.insert(59); OtherTFS.insert(65);
     OtherTFS.insert(68); OtherTFS.insert(72); OtherTFS.insert(78); OtherTFS.insert(79);
     // Look at NumEnergy randomly generated cases 
     srand(time(NULL));
     for (i=0; i<NumEnergy; i++)
     {
       cout << "State " << i << " ";
       for (j = 0; j < stretch_angle.size(); j++) {
	 if (OtherTFS.find(CapsomersNum[j]) != OtherTFS.end()) 
	 {
	   double randNumber = double(rand()%3);
	   stretch_angle[j]             = randNumber*M_PI/3.0;
	   VAngleShift[CapsomersNum[j]] = randNumber*M_PI/3.0; 
	   cout << randNumber << " ";
	 }
       }
       cout << endl;

       // assign new angle offset
       hex_body->setAngleOffset(stretch_angle);
	 
       // Solve 
       solver.solve(&model);
       double Wtot = bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy();
       cout << "State = " << i << " " << Wtot << endl;
       cout << "Wconf = " << hex_body->totalConformationalEnergy() << endl;
       
       // if ((i%100) == 0)
       // {
             SymmIndex = PrintTvirus.printMaster(i);
       // }

       cout << "Symmetry Index " << endl;
       cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
       cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
       cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
     }

   }
   else if(AnType == 2)
   {
     

   }






  



  


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
  cout << endl << "Center      = " << center  << endl;
  
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
  
  double  Wben = 0.0, Wpen = 0.0, Whex = 0.0, Wtot = 0.0, Wconf = 0.0;
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
  
  cout << " Max diff. in eta    = " << EtaMaxDiff       << endl;
  if (WcType == 5 || WcType == 6) {
    cout << " phiNode             = " << phiNode->point() << endl;
  }
  cout << " FvK number          = " << FvK              << endl;
  cout << " Bending energy      = " << Wben             << endl;
  cout << " Pent Stretch energy = " << Wpen             << endl;
  cout << " Hex  Stretch energy = " << Whex             << endl;
  cout << " Hex conform. energy = " << Wconf            << endl;
  cout << " Total energy        = " << Wtot             << endl;
  cout << " Total def. energy   = " << Wtot - Wconf     << endl;
  cout << " Asphericity         = " << asphericity      << endl;
  cout << " Ravg                = " << Ravg             << endl;
  cout << " AverageEdgeLength   = " << AvgEdgeLength    << endl;
  cout << "Final nodes symm index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
  cout << "Final eta   symm index (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
  cout << "Final theta symm index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
  cout << "Final nodes symm index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
  cout << "Final eta   symm index (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
  cout << "Final theta symm index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
  cout << "-----------------------------------------"   << endl;
  
  time (&end);
  dif = difftime (end,start);
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
