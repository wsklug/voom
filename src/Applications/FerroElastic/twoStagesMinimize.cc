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
#include "SkewMinHexonBody.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "MontecarloTwoStages.h"
#include "VoomMath.h"
//#include "print-body.cc"

using namespace voom;
using namespace std;

void printAllParaview(const std::string fileName, 
		      const Model::BodyContainer &, 
		      const vector<ScalarFieldNode<3>* > &, 
		      const vector<ScalarFieldNode<3>* > &,
		      const vector<unsigned int> &);

/*! This program creates a shell model with LoopShellBody for bending
  and 2 C0MembraneBodies for stretching, one of pentamers and one of
  hexamers.
 */
int main(int argc, char* argv[])
{
  time_t start,end;
  double dif;
  time (&start);

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
  // Parameters for MonteCarlo solver
  unsigned int maxIterMTS = 100;
  double finalRatioMTS = 1.0e-5;
  double T1 = 300.0;
  double T2 = 1.0;
  int VarEta = -2;
  int VarTheta = 0;
  
  string modelName;
  string parameterFileName = argv[1];
  ifstream inp, ifs;
 
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  inp >> modelName >> in_eta >> in_theta >> KC >> KG >> C0 >> mu >> kS >> min >> Ns >> maxIterMTS >> finalRatioMTS >> T1 >> T2 >> VarEta >> VarTheta;
  inp.close();

  // List input parameters
  cout << " modelName : " << modelName     << endl
       << " in_eta    : " << in_eta        << endl
       << " in_theta  : " << in_theta      << endl
       << " KC        : " << KC            << endl
       << " KG        : " << KG            << endl
       << " C0        : " << C0            << endl
       << " mu        : " << mu            << endl
       << " kS        : " << kS            << endl
       << " min       : " << min           << endl
       << " NshearDOF : " << Ns            << endl 
       << " MaxIter   : " << maxIterMTS    << endl
       << " FinalRatio: " << finalRatioMTS << endl
       << " T1        : " << T1            << endl
       << " T2        : " << T2            << endl
       << " VarEta    : " << VarEta        << endl
       << " VarTheta  : " << VarTheta      << endl; 
  
  
  
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



  // Read in the shear directions.  These are vectors in the current
  // configuration, pointing from vertex 1 to vertex 4 of each hexon,
  // as defined in Joe Rudnick's mathematica notebook. This vector
  // should point roughly along one of the edges of each triangular
  // element.
  //
  // Since we will define a "canonical" equilateral triangle reference
  // configuration for each element, we need to find which edge the
  // vector is closest to pointing along, and then choose that edge to
  // define the "shear angle" defining the reference pre-shear.
  
  while( token != "shear_direction") ifs >> token; 
  double AvgEdgeLength = 0.0, shear = 0.0, TOL = 2.5e-1, shear_dir_norm = 0.0;
  vector<double> shear_angle;
  vector< std::vector<int> > pent_connectivities, hex_connectivities;
  std::vector<int> cm(3,0.0);
  tvmet::Vector<double,3> shear_dir;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >> shear_dir(0) >> shear_dir(1) >> shear_dir(2);
    shear_dir_norm = tvmet::norm2(shear_dir); 
    
    for(j = 0; j < 3; j++) {
      cm[j] = connectivities[i](j);
    }
    
    // Edge vectors in current config.
    tvmet::Vector<double,3> e31(defNodes[cm[0]]->point()-defNodes[cm[2]]->point()), 
                            e32(defNodes[cm[1]]->point()-defNodes[cm[2]]->point()),
                            e12(defNodes[cm[1]]->point()-defNodes[cm[0]]->point());
    // Compute averate edge length for each triangle
    AvgEdgeLength += (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;

    // pentamers don't have shear directions
    if(shear_dir_norm == 0) 
    {
      pent_connectivities.push_back(cm);
    }
    // hexamers do
    else 
    {
      shear_dir = shear_dir/shear_dir_norm;
      hex_connectivities.push_back(cm);
      // Calculate the shear angle
      
      //normalize the edge vectors
      e31=e31/tvmet::norm2(e31);
      e32=e32/tvmet::norm2(e32);
      e12=e12/tvmet::norm2(e12);

      // Which edge is closest to the shear vector?  
      // Take cross product.  If zero (or < TOL) then two vectors are parallel or anti-parallel.  
      // Compute angles in degrees.

      if(tvmet::norm2(tvmet::cross(shear_dir,e31))<TOL) {

        if(tvmet::dot(shear_dir,e31)>0) shear=0.; // parallel with edge 31
        else shear=180.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(shear_dir,e32))<TOL){

        if(tvmet::dot(shear_dir,e32)>0) shear=60.; // parallel with 32
        else shear=240.; // anti-parallel

      } else if(tvmet::norm2(tvmet::cross(shear_dir,e12))<TOL) {

        if(tvmet::dot(shear_dir,e12)>0) shear=120.; // parallel with 12
        else shear=300.; // anti-parallel

      } 
      // If TOL is too small, it bombs.
      else {std::cout<<"Problem. The shear direction is not parallel to any triangle edge"<<std::endl; exit(1);}

      // Convert from degrees to radians.
      shear_angle.push_back(shear*PI/180.); 
    }
  }
  AvgEdgeLength /= ntri;
  cout << "Number of penton elements = " <<pent_connectivities.size() << endl
       << "Number of hexon elements = "  <<hex_connectivities.size() << endl;

  // Read capsomers to assign different shear dof to each one
  // Assume that hexons are written before pentons in the vtk input file.
  vector<unsigned int > CapsomersNum(ntri, 0.0);
  while( token != "capsomer") ifs >> token;
  ifs >> token;
  ifs >> token;
  ifs >> token;
  for (i = 0; i < ntri; i++)
  {
    ifs >> CapsomersNum[i];
    // cout << CapsomersNum[i] << endl;
  }
  cout << "Number of capsomers = " <<  CapsomersNum[ntri-1]+1 << endl;
  // Close mesh file in vtk format
  ifs.close();
  


  // Initialize shear node (ScalarFieldNode type)
  vector<ScalarFieldNode<3>* > shearNodes, directionNodes;
  shearNodes.reserve(Ns);   directionNodes.reserve(Ns);
  
  NodeBase::DofIndexMap idE(1);
  ScalarFieldNode<3>::PositionVector pE(0.0);
  for (i = 0; i < Ns; i++)
  {
    id++;
    idE[0] = dof++;
    ScalarFieldNode<3>* shearN = new ScalarFieldNode<3>(id, idE, pE, in_eta);
    shearNodes.push_back(shearN);
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
  LSB * bending_body = new LSB(bending, connectivities, nodes, quadOrder);
  // Use the body's output routines to generate vtk files with results
  bending_body->setOutput(paraview);
  bending_body->compute(true, false, false);
  cout << "En_Bend = " << bending_body->totalStrainEnergy() << endl;
  cout << "En_Bend = " << bending_body->energy() << endl;
  // bending_body->checkConsistency();

  

  // Create two stretching bodies, one for pents and one for hex's, which don't have any bending energy
  // (bending moduli set to zero)
  // Use EvansElastic Material, which separates area strain from pure shear strain (at constant area).  

  // Pentamers
  typedef EvansElastic PentStretchingMaterial;
  PentStretchingMaterial stretch(0.0, 0.0, 0.0, mu, kS);

  // Use generic C0MembraneBody for stretching of pents, with the EvansElastic material.
  typedef C0MembraneBody<TriangleQuadrature, PentStretchingMaterial, ShapeTri3> pent_membrane;
  pent_membrane * pentamer_body = new pent_membrane(stretch, pent_connectivities, nodes, quadOrder);
  pentamer_body->setOutput(paraview);
  pentamer_body->compute(true, false, false);
  cout << "En_Pent = " << pentamer_body->totalStrainEnergy() << endl;
  cout << "En_Pent = " << pentamer_body->energy() << endl;
  // pentamer_body->checkConsistency();
  
 
  
  // Hexamers
  // Use specialezed SkewHexonBody (modified version of C0MembraneBody), which has material hard-coded as
  // ModifiedEvansElasticMin, which takes care of the multiplicative F=AG "Eigen-strain" decomposition.

  // Correct the shear modulus if using double-well potential.
  // Not necessary if using single-well.
  // mu = mu/(2*(eta*eta+c2)+2*c1);
  TriangleQuadrature  HexQuad(1);
  const ShapeTri3::CoordinateArray dummy(0.0);
  ShapeTri3 HexShape(dummy);
  SkewMinHexonBody * hex_body = new SkewMinHexonBody(hex_connectivities, defNodes, shearNodes, directionNodes, shear_angle, mu, kS, HexQuad, &HexShape, CapsomersNum);
  hex_body->setOutput(paraview);
  hex_body->compute(true, false, false);
  cout << "En_Hex = " << hex_body->totalStrainEnergy() << endl;
  cout << "En_Hex = " << hex_body->energy() << endl;
  // hex_body->checkConsistency();
  

  
  // Set reference configuration
  // Define the edge length of the canonical 
  // const double EDGELEN = 2.0; - > AvgEdgeLength is used instead
  
  bending_body->SetRefConfiguration(AvgEdgeLength);
  pentamer_body->SetRefConfiguration(AvgEdgeLength);
  hex_body->SetRefConfiguration(AvgEdgeLength);
   
  //   pentamer_body->SetRefConfiguration(2.0);
  //   hex_body->SetRefConfiguration(2.0);


  
  

  // Create Model
  Model::BodyContainer bdc;
  bdc.push_back(bending_body);
  bdc.push_back(pentamer_body);
  bdc.push_back(hex_body);
  if (min == 1)   // shear magnitude is a dof, otherwise energy is minimized keeping eta fixed
  {
    nodes.insert(nodes.end(), shearNodes.begin(), shearNodes.end());
  }
  else if (min == 2)   // shear magnitude and direction are dof, otherwise energy is minimized keeping eta and theta fixed
  {
    nodes.insert(nodes.end(), shearNodes.begin(), shearNodes.end());
    nodes.insert(nodes.end(), directionNodes.begin(), directionNodes.end());
  }
  Model model(bdc, nodes);
  // Consistency check
  // model.checkConsistency(true,false);
  // model.checkRank(model.dof()-6,true);



  
  
  // Set solver and solve
  int iprint = 0;
  double factr = 1.0e1;
  double pgtol = 1.0e-3;
  int m = 5;
  int maxIter = 10000;
  // std::ifstream lbfgsbinp("lbfgsb.inp");
  // lbfgsbinp >> iprint >> factr >> pgtol >> m;
  cout << "Input iprint: " << iprint << endl
       << "Input factr:  " << factr  << endl
       << "Input pgtol:  " << pgtol  << endl
       << "Input m:      " << m      << endl;
  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);
  // Solve to minimize the energy
  // bdc[0]->printParaview("initial");
  printAllParaview("Initial.vtk", bdc, shearNodes, directionNodes, CapsomersNum);
  // solver.solve(&model);

  vector<int > VarType;
  vector<vector<ScalarFieldNode<3>* > > MontecarloDof;
  if (VarTheta > -2)
  {
    MontecarloDof.push_back(directionNodes);
    VarType.push_back(VarTheta);
  }
  if (VarEta > -2)
  {
    MontecarloDof.push_back(shearNodes);
    VarType.push_back(VarEta);
  }
  MontecarloTwoStages MTS(MontecarloDof, VarType, &solver, maxIterMTS, true); 
  MTS.SetTempSchedule(MTS.EXPONENTIAL, T1, T2, finalRatioMTS);
  MTS.solve(&model); 

  printAllParaview("Minimized.vtk", bdc, shearNodes, directionNodes, CapsomersNum);

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


  // Compute max difference in eta
  double EtaMaxDiff = 0.0, EtaTempDiff = 0.0;
  if (Ns > 1)
  {
    EtaMaxDiff = abs(shearNodes[1]->point()-shearNodes[0]->point());
    for (i = 1; i < Ns-1; i++)
    {
      EtaTempDiff = abs(shearNodes[i+1]->point()-shearNodes[i]->point());
      if ( EtaTempDiff > EtaMaxDiff ) EtaMaxDiff = EtaTempDiff;
    }
  }

  //---------------//
  // Print results //
  // Print the total energies of three bodies
  
  double FvK = 0.0, Y = 0.0, Wben = 0.0, Wpen = 0.0, Whex = 0.0, Wtot = 0.0;
  Y = 4.0*kS*mu/(kS+mu);
  FvK = Y*Ravg*Ravg/KC;
  Wben =  bending_body->totalStrainEnergy();
  Wpen =  pentamer_body->totalStrainEnergy();
  Whex =  hex_body->totalStrainEnergy();
  Wtot = Wben + Wpen + Whex;
  cout << "----------------------------------------" << endl;
  if (shearNodes.size() == 1)
  {
    cout << " Shear magnitude     = " << shearNodes[0]->point()    << endl;
    cout << " Direction angle     = " << directionNodes[0]->point()<< endl;
  }
  else
  {
    cout << " Shear magnitude     = " << endl;
    for (i = 0; i < Ns; i++)
    {
      cout << shearNodes[i]->point() << "  ";
      if ( (i+1)%5 == 0 || i == Ns-1 ) cout << endl;
    }

    cout << " Direction angle     = " << endl;
    for (i = 0; i < Ns; i++)
    {
      cout << directionNodes[i]->point() << "  ";
      if ( (i+1)%5 == 0 || i == Ns-1 ) cout << endl;
    }
  }
  cout << " Max diff. in eta    = " << EtaMaxDiff    << endl;
  cout << " FvK number          = " << FvK           << endl;
  cout << " Bending energy      = " << Wben          << endl;
  cout << " Pent Stretch energy = " << Wpen          << endl;
  cout << " Hex  Stretch energy = " << Whex          << endl;
  cout << " Total energy        = " << Wtot          << endl;
  cout << " Asphericity         = " << asphericity   << endl;
  cout << " Ravg                = " << Ravg          << endl;
  cout << " AverageEdgeLength   = " << AvgEdgeLength << endl;
  cout << "----------------------------------------" << endl;
  
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

  if (min==0)
  {
    for (i = 0; i<Ns; i++)
    {
      delete shearNodes[i];
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


void printAllParaview(const std::string fileName,
		      const Model::BodyContainer & bdc, 
		      const vector<ScalarFieldNode<3>* > & shearNodes, 
		      const vector<ScalarFieldNode<3>* > & directionNodes,  
		      const vector<unsigned int> & CapsomersNum)
{
  // OutPut bending and stretching energy in current configuration
  ofstream ofs(fileName.c_str());
  if (!ofs) {
    std::cout << "can not open output ("
	      << fileName
	      << ") file." << std::endl;
    exit(0);
  }

  // Start from bending body which has all elements
  Body::ElementContainer BendingElements = bdc[0]->elements();
  Body::NodeContainer BendingNodes = bdc[0]->nodes();

  // Node Section
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << endl
	<< "ASCII" << endl
	<< "DATASET POLYDATA" << endl
        << "POINTS  " << BendingNodes.size() << "  double" << endl;

  // Output nodal postions
  Body::NodeIterator pn; 
  for (pn = BendingNodes.begin(); pn!= BendingNodes.end(); pn ++)
  {
    DeformationNode<3> * node = dynamic_cast<DeformationNode<3>* >(*pn);

    if (node == NULL) exit(0);

    const Vector3D & nodalPos =  node->point();
    ofs << std::setprecision(16) 
	<< nodalPos(0) << "  "
	<< nodalPos(1) << "  "
	<< nodalPos(2) << std::endl;
  }

  // Element Section
  int nel = BendingElements.size();
  ofs << "POLYGONS  " << nel << "  "
      << 4*nel << endl;
  for(int e = 0; e < nel; e++)
  {
    if(!(bdc[0]->active(e)) ) exit(0); // All elements should be active in this application
    // ::FElement_t * pe = dynamic_cast<LoopShellBody::FElement_t *>(BendingElements[e]);
    ofs << 3 << "  "
	<< setw(10) << BendingElements[e]->baseNodes()[0] -> id()
	<< setw(10) << BendingElements[e]->baseNodes()[1] -> id()
	<< setw(10) << BendingElements[e]->baseNodes()[2] -> id()
	<< endl;
  }

  // Cell section
  ofs << "CELL_DATA    " << nel << endl;
  // Bending energy
  ofs << "SCALARS    BendingEnergy    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  for(int e = 0; e < nel; e++) {
      ofs << BendingElements[e]->energy() << endl;
  }
  ofs << std::endl;

  // Stretching energy
  Body::ElementContainer PentElements = bdc[1]->elements();
  Body::ElementContainer HexElements = bdc[2]->elements();
  ofs << "SCALARS    StretchingEnergy    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  for(int e = 0; e < HexElements.size(); e++) {
      ofs << HexElements[e]->energy() << endl;
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << PentElements[e]->energy() << endl;
  }
  ofs << endl;

  // Shear magnitude
  ofs << "SCALARS    ShearMagnitude    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  if (shearNodes.size() == 1) {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << shearNodes[0]->point() << endl;
    }
  }
  else {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << shearNodes[CapsomersNum[e]]->point() << endl;
    }
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;

  // Shear direction
  ofs << "SCALARS    ShearDirection    double    1" << endl;
  ofs << "LOOKUP_TABLE default" << endl;
  
  if (directionNodes.size() == 1) {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << directionNodes[0]->point() << endl;
    }
  }
  else {
    for(int e = 0; e < HexElements.size(); e++) {
      ofs << directionNodes[CapsomersNum[e]]->point() << endl;
    }
  }

  for(int e = 0; e < PentElements.size(); e++) {
      ofs << 0.0 << endl;
  }
  ofs << endl;

  ofs.close();

}
