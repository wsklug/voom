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
#include "SkewHexonBody.h"
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
                          // min = 4 Special case when phi is not a DOF with WcType = 5 or 6, theta is DOF, eta is not a DOF
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
  // Additional energy input parameters
  int WcType = 1;
  double WcConst0 = 1.0;
  double WcConst1 = 1.0;
  double WcConst2 = 1.0;
  double WcConst3 = 1.0;
  double WcConst[4];
  int WcStepsEta = 1;
  int WcStepsTheta = 1;

  int Tnumber = 0;

  // Output
  unsigned int NPout = 0; // Nodes needed to build the outline
  double HexonShift = 0.0;// To separate direction arrows
  string OutFile;
  int refinement = 1;

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
  in_theta *= (M_PI/180.0);   // in_theta input is in degrees, need to be transformed in radiants
  inp >> dump >> in_phi;
  in_phi *= (M_PI/180.0);     // in_phi input is in degrees, need to be transformed in radiants
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
  inp >> dump >> WcStepsEta;
  inp >> dump >> WcStepsTheta;
  inp >> dump >> NPout;
  inp >> dump >> HexonShift;
  inp >> dump >> OutFile;
  inp >> dump >> refinement;
  inp >> dump >> Tnumber;
  
  inp.close();

  // Populate Wconst
  WcConst[0] = WcConst0;
  WcConst[1] = WcConst1;
  WcConst[2] = WcConst2;
  WcConst[3] = WcConst3;

  // List input parameters
  cout << " modelName    : " << modelName    << endl
       << " OutlineName  : " << OutlineName  << endl
       << " min          : " << min          << endl
       << " Ns_eta       : " << Ns_eta       << endl
       << " Ns_theta     : " << Ns_theta     << endl
       << " in_eta       : " << in_eta       << endl
       << " in_theta     : " << in_theta     << endl
       << " in_phi       : " << in_phi       << endl
       << " KC           : " << KC           << endl
       << " KG           : " << KG           << endl
       << " C0           : " << C0           << endl
       << " mu           : " << mu           << endl
       << " kS           : " << kS           << endl
       << " WcType       : " << WcType       << endl
       << " alpha        : " << WcConst[0]   << endl
       << " beta         : " << WcConst[1]   << endl
       << " gamma        : " << WcConst[2]   << endl 
       << " delta        : " << WcConst[3]   << endl 
       << " WcStepsEta   : " << WcStepsEta   << endl
       << " WcStepsTheta : " << WcStepsTheta << endl
       << " NPout        : " << NPout        << endl
       << " HexonShift   : " << HexonShift   << endl
       << " OutFile      : " << OutFile      << endl
       << " Refinement   : " << refinement   << endl
       << " Tnumber      : " << Tnumber      << endl;
      
  



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
  // 4 -> theta but phi is not a DOF
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
    case 4:
      NumDOF = nodes.size()*3 + Ns_theta;
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



  // For T4
     vector<double > VAngleShift(directionNodes.size(), 0.0);
     if (Tnumber == 4)
     {
       VAngleShift[0]  = M_PI/3.0;  VAngleShift[1]  =-M_PI/3.0;  VAngleShift[2] = 0.0;        VAngleShift[3] =-M_PI/3.0;  VAngleShift[4] = 0.0;
       VAngleShift[5]  = M_PI/3.0;  VAngleShift[6]  =-M_PI/3.0;  VAngleShift[7] = M_PI/3.0;   VAngleShift[8] =-M_PI/3.0;  VAngleShift[9] = M_PI/3.0;
       VAngleShift[10] =-M_PI/3.0;  VAngleShift[11] = 0.0;       VAngleShift[12] =-M_PI/3.0;  VAngleShift[13] = 0.0;      VAngleShift[14] = 0.0;
       
       VAngleShift[15] =-M_PI/3.0;  VAngleShift[16] = 0.0;       VAngleShift[17] = 0.0;       VAngleShift[18] = M_PI/3.0; VAngleShift[19] =-M_PI/3.0;
       VAngleShift[20] =-M_PI/3.0;  VAngleShift[21] = 0.0;       VAngleShift[22] = 0.0;       VAngleShift[23] = 0.0;      VAngleShift[24] = M_PI/3.0;
       VAngleShift[25] =-M_PI/3.0;  VAngleShift[26] =-M_PI/3.0;  VAngleShift[27] = M_PI/3.0;  VAngleShift[28] = 0.0;      VAngleShift[29] = M_PI/3.0;
  

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
     }

  // For T9
     if (Tnumber == 9)
     {
       // Set the initial pre-stretch configuration in the same configuration as the minimum 
       // configuration found using InfCur
       // Shift in initial hexamer configuration
       inp.open("T9_23Interactions.dat", ios::in);
       if (!inp) {
	 cout << "Cannot open input file: T9_23Interactions.dat" << endl;
	 return(0);
       }
  
       for (i = 0; i < VAngleShift.size(); i++)
       {
	 inp >> VAngleShift[i];
	 cout << VAngleShift[i] << " ";
	 if ( (i+1)%5 == 0 || i == VAngleShift.size()-1 ) cout << endl;
	 VAngleShift[i] *= (M_PI/3.0);
       }

       unsigned int CurrentCap = 0;
       uint ind = 0;
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
       // End of setting initial configuration



       // Randomize three fold sites in initial configuration
       // Other Three fold site hexamers
       set<uint > OtherTFS;
       OtherTFS.insert(3);  OtherTFS.insert(8);  OtherTFS.insert(24); OtherTFS.insert(73);
       OtherTFS.insert(15); OtherTFS.insert(20); OtherTFS.insert(26); OtherTFS.insert(31);
       OtherTFS.insert(35); OtherTFS.insert(39); OtherTFS.insert(43); OtherTFS.insert(49);
       OtherTFS.insert(52); OtherTFS.insert(56); OtherTFS.insert(59); OtherTFS.insert(65);
       OtherTFS.insert(68); OtherTFS.insert(72); OtherTFS.insert(78); OtherTFS.insert(79);

       // Select 20 random numbers (equivalent to the number of other 3 fold sites)
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
       // End of randomize backgroud
  
       
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
  case 3: // stretch direction and magnitude are a dof
    nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
    if (WcType == 5 || WcType == 6) {
      nodes.push_back(phiNode);
    }
    break;
  case 4: // stretch direction is a dof (but phi is not)
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
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
  int m = 5;             // Previously 10
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


   // Solve to minimize the energy
   if (WcType == 0)
   {
     // First solve with No conformational energy
     solver.solve(&model);
     cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
     cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
       cout << " Stretch magnitude     = " << endl;
       for (j = 0; j < Ns_eta; j++)
       {
	 cout << stretchNodes[j]->point() << "  ";
	 if ( (j+1)%5 == 0 || j == Ns_eta-1 ) cout << endl;
       }
       cout << " Direction angle     = " << endl;
       for (j = 0; j < Ns_theta; j++)
       {
	 cout << directionNodes[j]->point() << "  ";
	 if ( (j+1)%5 == 0 || j == Ns_theta-1 ) cout << endl;
       }
       SymmIndex = PrintTvirus.printMaster(1);
       cout << "Symmetry Index " << endl;
       cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
       cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
       cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
       cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
       cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
       cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;


       if (Tnumber == 9)
       {
	 double randNumber = double(rand()%2)+1.0;
	 cout << "randNumber = " << randNumber << endl;
	 VAngleShift[26] += randNumber*M_PI/3.0; 
	 VAngleShift[43] += randNumber*M_PI/3.0; 
	 VAngleShift[59] += randNumber*M_PI/3.0; 
	 for (j = 0; j < stretch_angle.size(); j++) {
	   if (CapsomersNum[j] == 26) 
	   {
	     stretch_angle[j]             += randNumber*M_PI/3.0;
	   }
	   else if (CapsomersNum[j] == 43) 
	   {
	     stretch_angle[j]             += randNumber*M_PI/3.0; 
	   }
	   else if (CapsomersNum[j] == 59) 
	   {
	     stretch_angle[j]             += randNumber*M_PI/3.0; 
	   }
	 }
	 // assign new angle offset
	 hex_body->setAngleOffset(stretch_angle);
	 // Solve 
	 solver.solve(&model);
	 double Wtot = bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy();
	 cout<< "Wtot = " << Wtot << endl;

	 SymmIndex = PrintTvirus.printMaster(2);
	 cout << "Symmetry Index " << endl;
	 cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	 cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	 cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	 cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	 cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	 cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
       }
       
   }



   if (WcType == 3 || WcType == 5 || WcType == 6)
   {
     // First solve with No conformational energy
     WcConst[0] = 0.0;
     solver.solve(&model);
         cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
         cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
	 cout << " Stretch magnitude     = " << endl;
	 for (j = 0; j < Ns_eta; j++)
	 {
	   cout << stretchNodes[j]->point() << "  ";
	   if ( (j+1)%5 == 0 || j == Ns_eta-1 ) cout << endl;
	 }
	 cout << " Direction angle     = " << endl;
	 for (j = 0; j < Ns_theta; j++)
	 {
	   cout << directionNodes[j]->point() << "  ";
	   if ( (j+1)%5 == 0 || j == Ns_theta-1 ) cout << endl;
	 }
	 if (WcType == 5 || WcType == 6) {
	   cout << " phiNode             = " << phiNode->point() << endl;
	 }
	 SymmIndex = PrintTvirus.printMaster(1);
	 cout << "Iter 1 - Symmetry indices " << endl;
	 cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	 cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	 cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	 cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	 cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	 cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;

     // Then apply conformational energy to apply bias in theta
     WcConst[0] = WcConst0;
     for (i = 1; i < WcStepsTheta+1; i++)
     {
       solver.solve(&model);
           cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
	   cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
	   cout << " Stretch magnitude     = " << endl;
	   for (j = 0; j < Ns_eta; j++)
	   {
	     cout << stretchNodes[j]->point() << "  ";
	     if ( (j+1)%5 == 0 || j == Ns_eta-1 ) cout << endl;
	   }
	   cout << " Direction angle     = " << endl;
	   for (j = 0; j < Ns_theta; j++)
	   {
	     cout << directionNodes[j]->point() << "  ";
	     if ( (j+1)%5 == 0 || j == Ns_theta-1 ) cout << endl;
	   }
	   if (WcType == 5 || WcType == 6) {
	     cout << " phiNode             = " << phiNode->point() << endl;
	   }
	   SymmIndex = PrintTvirus.printMaster(i+1);
	   cout << "Iter " << i+1 << " - Symmetry indices " << endl;
           cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	   cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	   cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	   cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	   cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	   cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
	   cout << "WcConst = [ " << WcConst[0] << " , " << WcConst[1] << " ] " << endl;
	   WcConst[0] *= 2.0; // i*WcConst0; 
     }

     if (Tnumber == 4)
     {
       // Then apply gradually the bias in eta
       for (i = 1; i < WcStepsEta+1; i++)
       {
	 in_eta += 0.01;
	 for (j = 0; j < Ns_eta; j++) {
	   stretchNodes[j]->setPoint(in_eta);
	 }
	 
	 solver.solve(&model);
	     cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
	     cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
	     cout << " Stretch magnitude     = " << endl;
	     for (j = 0; j < Ns_eta; j++)
	     {
		 cout << stretchNodes[j]->point() << "  ";
		 if ( (j+1)%5 == 0 || j == Ns_eta-1 ) cout << endl;
	     }
	     cout << " Direction angle     = " << endl;
	     for (j = 0; j < Ns_theta; j++)
	     {
		 cout << directionNodes[j]->point() << "  ";
		 if ( (j+1)%5 == 0 || j == Ns_theta-1 ) cout << endl;
	     }
	     if (WcType == 5 || WcType == 6) {
	       cout << " phiNode             = " << phiNode->point() << endl;
	     }
	     SymmIndex = PrintTvirus.printMaster(i+1+WcStepsTheta);
	     cout << "Iter " << i+1 << " - Symmetry indices " << endl;
	     cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	     cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	     cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	     cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	     cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	     cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
       }
     }



   }

   




   if (WcType == 4)
   {
     // First bias eta toward value set in WcConst[3]
     WcConst[0] = 0.0;
     for (i = 1; i < WcStepsEta+1; i++)
     {
       WcConst[2] = i*WcConst2; 
       solver.solve(&model);
       cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
       cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
       SymmIndex = PrintTvirus.printMaster(i);
       cout << "Iter " << i+1 << " - Symmetry indices " << endl;
           cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	   cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	   cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	   cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	   cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	   cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;

       cout << " Stretch magnitude     = " << endl;
       for (j = 0; j < Ns_eta; j++)
       {
	 cout << stretchNodes[j]->point() << "  ";
	 if ( (j+1)%5 == 0 || j == Ns_eta-1 ) cout << endl;
       }
     }
   
     // Second force positive stretch to align vertex to vertex
     for (i = 1; i < WcStepsTheta+1; i++)
     {
       WcConst[0] = i*WcConst0; 
       solver.solve(&model);
       cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
       cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
       SymmIndex = PrintTvirus.printMaster(i+WcStepsEta);
       cout << "Iter " << i+1 << " - Symmetry indices " << endl;
           cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	   cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	   cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	   cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	   cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	   cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
     }

     // Third remove bias to eta
     for (i = 1; i < WcStepsEta+1; i++)
     {
       WcConst[2] -= WcConst2; 
       solver.solve(&model);
       cout << "Wtot  = " << bending_body->totalStrainEnergy() +  pentamer_body->totalStrainEnergy() + hex_body->totalStrainEnergy() << endl;
       cout << "Wconf = " <<  hex_body->totalConformationalEnergy() << endl;
       SymmIndex = PrintTvirus.printMaster(i+WcStepsEta+WcStepsTheta);
       cout << "Iter " << i+1 << " - Symmetry indices " << endl;
           cout << "Nodes index (hex) = " << std::setprecision(12) << SymmIndex[0] << endl;
	   cout << "Eta index   (hex) = " << std::setprecision(12) << SymmIndex[1] << endl;
	   cout << "Theta index (hex) = " << std::setprecision(12) << SymmIndex[2] << endl;
	   cout << "Nodes index (el)  = " << std::setprecision(12) << SymmIndex[3] << endl;
	   cout << "Eta index   (el)  = " << std::setprecision(12) << SymmIndex[4] << endl;
	   cout << "Theta index (el)  = " << std::setprecision(12) << SymmIndex[5] << endl;
     }
     cout << "WcConst = [ " << WcConst[0] << " , " << WcConst[1] << " , " << WcConst[2] << " , " << WcConst[3] <<" ] " << endl;
   }


   


  // Printing results
  // Print the postprocessing strains and stress
  //hex_body->printParaviewEigVec("Hexbody");
  //pentamer_body->printParaviewEigVec("Pentbody");
  //hex_body->printParaview2("Hexbody2");
  //pentamer_body->printParaview2("Pentbody2");





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
