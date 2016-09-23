#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
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
#include "VoomMath.h"

using namespace voom;
using namespace std;

#include "UtilsConformationalChanges.h"



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
  
  // string InputFileName = argv[1];
  
 
  
  

  InputData InpVal ReadInputFile(argv[1]);

  

  // Read in a series of mesh files with different shearing directions
  // Loop over possible intial direction configurations
  ofstream ResSum("ResultsSummary.dat");
  if (!ResSum) {
    std::cout << "can not open output ResultsSummary.dat" << std::endl;
    exit(0);
  }

  ifstream ifs;
  for (unsigned int symmConf = 1; symmConf < 28; symmConf++)
  {
    // Read in the mesh, in (ascii) legacy vtk format.
    stringstream MeshFile;
    string meshFileName;
    MeshFile << "/u/home/campus/luigiemp/ConformationalChange/T16mesh/T16vec3s" << symmConf << ".vtk";
    meshFileName = MeshFile.str();

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



    // Read in the stretch directions.  These are vectors in the current
    // configuration, pointing from vertex 1 to vertex 4 of each hexon,
    // as defined in Joe Rudnick's mathematica notebook. This vector
    // should point roughly along one of the edges of each triangular
    // element.
    //
    // Since we will define a "canonical" equilateral triangle reference
    // configuration for each element, we need to find which edge the
    // vector is closest to pointing along, and then choose that edge to
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
	    else {std::cout<<"Problem. The stretch direction is not parallel to any triangle edge"<<std::endl; exit(1);}

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
	    stretch_angle.push_back((stretch+AngleShift)*PI/180.); 
	  }
      }
    AvgEdgeLength /= ntri;
    cout << "Number of penton elements = " <<pent_connectivities.size() << endl
	 << "Number of hexon elements = "  <<hex_connectivities.size() << endl;

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
	// if (CapsomersNum[i] > 95 && CapsomersNum[i] < 100) cout << CapsomersNum[i] << endl;
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

    // Correct the stretch modulus if using double-well potential.
    // Not necessary if using single-well.
    // mu = mu/(2*(eta*eta+c2)+2*c1);
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
    if (min == 1)   // stretch magnitude is a dof, otherwise energy is minimized keeping eta fixed
      {
	nodes.insert(nodes.end(), stretchNodes.begin(), stretchNodes.end());
      }
    else if (min == 2)   // stretch magnitude and direction are dof, otherwise energy is minimized keeping eta and theta fixed
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
    std::ifstream lbfgsbinp("lbfgsb.inp");
    lbfgsbinp >> iprint >> factr >> pgtol >> m;
    cout << "Input iprint: " << iprint << endl
	 << "Input factr:  " << factr  << endl
	 << "Input pgtol:  " << pgtol  << endl
	 << "Input m:      " << m      << endl;
  
    // Solve to minimize the energy
    Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter);
  
      
      
      stringstream OutA, OutB, OutC, OutD, OutK;
      string OutFileName, FileNameB;
      // Print initial results to vtk files
      OutA << "Initial_" << symmConf << ".vtk";
      OutFileName = OutA.str();
      OutK << "T16outline_" << symmConf << ".vtk";
      FileNameB = OutK.str();
      printAllParaview(OutFileName, FileNameB, OutlineName, NPout, bdc, stretchNodes, directionNodes, CapsomersNum);

      OutB << "StretchInitialDirections_" << symmConf << ".vtk";
      OutFileName = OutB.str();
      printStretchDirectioParaview(OutFileName, meshFileName, bdc, AngleShift, HexonShift, stretchNodes, directionNodes, CapsomersNum, false);

      // Solve
      solver.solve(&model);


      OutC << "Minimized_" << symmConf << ".vtk";
      OutFileName = OutC.str();
      printAllParaview(OutFileName, FileNameB, OutlineName, NPout, bdc, stretchNodes, directionNodes, CapsomersNum);
     
      OutD << "StretchFinalDirections_" << symmConf << ".vtk";
      OutFileName = OutD.str();
      printStretchDirectioParaview(OutFileName, meshFileName, bdc, AngleShift, HexonShift, stretchNodes, directionNodes, CapsomersNum, true);


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
      ResSum << "----------------------------------------  " << symmConf << endl;
      ResSum << "Center      = " << center  << endl;
  

      for(i = 0; i < defNodes.size(); i++)
	{
	  Rtemp = tvmet::norm2(defNodes[i]->point());
	  deltaR += pow(Rtemp-Ravg,2);
	}
      asphericity = sqrt( deltaR/defNodes.size() )/Ravg;


      // Compute max difference in eta and theta
      double EtaMaxDiff = 0.0, EtaTempDiff = 0.0, ThetaMaxDiff = 0.0, ThetaTempDiff = 0.0;
      if (Ns > 1)
	{
	  EtaMaxDiff = abs(stretchNodes[1]->point()-stretchNodes[0]->point());
	  ThetaMaxDiff = abs(directionNodes[1]->point()-directionNodes[0]->point());
	  for (i = 1; i < Ns-1; i++)
	    {
	      EtaTempDiff = abs(stretchNodes[i+1]->point()-stretchNodes[i]->point());
	      ThetaTempDiff = abs(directionNodes[i+1]->point()-directionNodes[i]->point());
	      if ( EtaTempDiff > EtaMaxDiff ) EtaMaxDiff = EtaTempDiff;
	      if ( ThetaTempDiff > ThetaMaxDiff ) ThetaMaxDiff = ThetaTempDiff;
	    }
	}

      //---------------//
      // Print results //
      // Print the total energies of three bodies
  
      double FvK = 0.0, Y = 0.0, Wben = 0.0, Wpen = 0.0, Whex = 0.0, Wtot = 0.0, Wconf = 0.0;
      Y = 4.0*kS*mu/(kS+mu);
      FvK = Y*Ravg*Ravg/KC;
      Wben  =  bending_body->totalStrainEnergy();
      Wpen  =  pentamer_body->totalStrainEnergy();
      Whex  =  hex_body->totalStrainEnergy();
      Wconf =  hex_body->totalConformationalEnergy();
      Wtot = Wben + Wpen + Whex;
      cout << "----------------------------------------" << endl;
      if (stretchNodes.size() == 1)
	{
	  cout << " Stretch magnitude     = " << stretchNodes[0]->point()    << endl;
	  ResSum << " Stretch magnitude     = " << stretchNodes[0]->point()    << endl;
	  cout << " Direction angle     = " << directionNodes[0]->point()<< endl;
	  ResSum << " Direction angle     = " << directionNodes[0]->point()<< endl;
	}
      else
	{
	  cout << " Stretch magnitude     = " << endl;
	  for (i = 0; i < Ns; i++)
	    {
	      cout << stretchNodes[i]->point() << "  ";
	      ResSum << stretchNodes[i]->point() << "  ";
	      if ( (i+1)%5 == 0 || i == Ns-1 ) {
		cout << endl;
		ResSum << endl;
	      }
	    }
	  
	  cout << " Direction angle     = " << endl;
	  ResSum << " Direction angle     = " << endl;
	  for (i = 0; i < Ns; i++)
	    {
	      cout << directionNodes[i]->point() << "  ";
	      ResSum << directionNodes[i]->point() << "  ";
	      if ( (i+1)%5 == 0 || i == Ns-1 ) {
		cout << endl;
		ResSum << endl;
	      }
	    }
	}
      cout << " Max diff. in eta    = " << EtaMaxDiff    << endl;
      cout << " Max diff. in theta  = " << ThetaMaxDiff  << endl;
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

      ResSum << " Max diff. in eta    = " << EtaMaxDiff    << endl;
      ResSum << " Max diff. in theta  = " << ThetaMaxDiff  << endl;
      ResSum << " FvK number          = " << FvK           << endl;
      ResSum << " Bending energy      = " << Wben          << endl;
      ResSum << " Pent Stretch energy = " << Wpen          << endl;
      ResSum << " Hex  Stretch energy = " << Whex          << endl;
      ResSum << " Hex conform. energy = " << Wconf         << endl;
      ResSum << " Total energy        = " << Wtot          << endl;
      ResSum << " Total def. energy   = " << Wtot - Wconf  << endl;
      ResSum << " Asphericity         = " << asphericity   << endl;
      ResSum << " Ravg                = " << Ravg          << endl;
      ResSum << " AverageEdgeLength   = " << AvgEdgeLength << endl;
      ResSum << "----------------------------------------" << endl << endl << endl << endl;

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
          
  } // end of loop over different directions
  ResSum.close();



  time (&end);
  dif = difftime (end,start);
  cout << endl << "All done :) in " << dif  << " s" << endl;

  return (0);
  
}


