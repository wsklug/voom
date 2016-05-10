#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "CGfast.h"
#include "Contact.h"
#include "ViscousRegularizer.h"
#include "RigidHemisphereAL.h"
#include "RigidPlateAL.h"

#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCell.h>

#include "Morse.h"
#include "SpringPotential.h"
#include "PotentialBody.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void insertValenceInVtk(const std::string fileName,\
			vtkSmartPointer<vtkPolyData> mesh);

std::vector<double> calcEdgeLenAndStdDev
(std::vector< DeformationNode<3>* > a, 
 vector< tvmet::Vector<int,3> > b);

int main(int argc, char* argv[])
{
  clock_t t1,t2;
  t1=clock();
  if( argc < 2 ) {
    cout << "Usage: indent modelName [-g gamma -p prestressFlag -n nsteps]."
	 << endl;
    return(0);
  }
  
#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif
  
  bool verbose=true;
#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif
  
  for(int i=0; i<argc; i++) {
    std::cout << std::setw(8) << i << "\t"
	      << argv[i] << std::endl;
  }

  string modelName = argv[1];

  double gamma_inp=0.0;
  double indent_inp=0.0;
  double viscosity_inp=0.0;
  double friction_inp=0.0;
  double minStep_inp=0.01;
  double step_inp=0.0;
  bool remesh = true;
  string prestressFlag = "yes";
  bool badCommandLine=false;
  bool unload=false;
  int option;

  optind=2; // start at argv[2] looking for command-line options
  while( (option=getopt(argc,argv,"e:f:g:i:m:p:s:v:u:")) != -1 ) {
    std::cout << "option = " << char(option) << std::endl;
    switch (option) {
    case 'g':
      gamma_inp = std::atof(optarg);
      std::cout << "gamma: " << gamma_inp << std::endl;
      break;
    case 'i':
      indent_inp = std::atof(optarg);
      std::cout << "max indentation: " << indent_inp << std::endl;
      break;
    case 'm':
      minStep_inp = std::atof(optarg);
      std::cout << "min step: " << minStep_inp << std::endl;
      break;
    case 'p' : 
      prestressFlag = std::string(optarg);
      std::cout << "prestressFlag: " << prestressFlag << std::endl;
      break;
    case 's':
      step_inp = std::atof(optarg);
      std::cout << "step size: " << step_inp << std::endl;
      break;
    case 'v':
      viscosity_inp = std::atof(optarg);
      std::cout << "viscosity: " << viscosity_inp << std::endl;
      break;
    case 'f':
      friction_inp = std::atof(optarg);
      std::cout << "friction: " << friction_inp << std::endl;
      break;
    case 'u':
      unload = true;
      std::cout << "Unloading curve will be calculated." << std::endl;
      break;
    case '?':
      std::cout << "Option " << char(optopt) << " is not supported." 
		<< std::endl;
      return 1;
    default :
      badCommandLine=true;
      break;
    }
  }

  std::cout << "optind = " << optind << std::endl
	    << "option = "<< option << std::endl;

  if(gamma_inp <=0.0) {
    std::cout << "gamma = " << gamma_inp << " but should be positive." 
	      << std::endl;
    return 0;
  }

  //For Morse material 
  double epsilon;
  double percentStrain;
  double pressureFactor;
  bool harmonicRelaxNeeded;
  int interimIter = 10;
  int continueFromNum = 1;
  double Ravg = 0;
  double Rshift = 1.0;
  double Zmin;
  double Zmax;
  double afmR = 1.0*Ravg;
  tvmet::Vector<double,3> xc(0.0);
  double Z_glass;
  double dZ;

  //Read epsilon and percentStrain from input file. percentStrain is
  //calculated so as to set the inflection point of Morse potential
  //at a fixed distance relative to the equilibrium separation
  //e.g. 1.1*R_eq, 1.5*R_eq etc.
  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  string temp;
  miscInpFile >> temp >> epsilon
	      >> temp >> percentStrain
	      >> temp >> pressureFactor
	      >> temp >> harmonicRelaxNeeded
	      >> temp >> interimIter
	      >> temp >> continueFromNum
	      >> temp >> Ravg
	      >> temp >> Rshift
	      >> temp >> Zmax
	      >> temp >> Zmin
	      >> temp >> Z_glass
	      >> temp >> dZ
	      >> temp >> xc[0]
	      >> temp >> xc[1]
	      >> temp >> xc[2];
  
  miscInpFile.close();

  vtkDataSetReader * reader = vtkDataSetReader::New();
  vtkSmartPointer<vtkPolyData> mesh;
  string inputFileName;
  vtkSmartPointer<vtkDataArray> displacements;

  vtkSmartPointer<vtkPolyData> mesh_prev;
  vtkSmartPointer<vtkDataArray> displacements_prev;

  if(continueFromNum > 0){
    harmonicRelaxNeeded = false;
    char name[100]; 
    
    //Read the latest VTK file
    sprintf(name,"%s-body%d-step%04d",modelName.c_str(),1,continueFromNum);
    inputFileName = string(name);
    reader->SetFileName( inputFileName.c_str() );
    mesh = reader->GetPolyDataOutput();
    mesh->Update();
    displacements = mesh->GetPointData()->GetVectors("displacements");

    //Also read the second to last VTK file
    sprintf(name,"%s-body%d-step%04d",modelName.c_str(),1,continueFromNum-1);
    inputFileName = string(name);
    reader->SetFileName( inputFileName.c_str() );
    mesh_prev = reader->GetPolyDataOutput();
    mesh_prev->Update();
    displacements_prev = mesh->GetPointData()->GetVectors("displacements");
    

  }
  else{
    inputFileName = modelName + ".vtk";
    reader->SetFileName( inputFileName.c_str() );
    //We will use this object, shortly, to ensure consistent triangle orientations
    vtkPolyDataNormals * normals = vtkPolyDataNormals::New();

    //We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput()
    //If our input vtk file has vtkUnstructuredGridData instead of vtkPolyData
    //then we need to convert it using vtkGeometryFilter
    vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();
    ds->Update();
    if(ds->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
      vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = 
	reader->GetUnstructuredGridOutput();    
      vtkSmartPointer<vtkGeometryFilter> geometryFilter = 
	vtkSmartPointer<vtkGeometryFilter>::New();
      geometryFilter->SetInput(unstructuredGrid);
      geometryFilter->Update(); 
      vtkSmartPointer<vtkPolyData> polydata = geometryFilter->GetOutput();
      normals->SetInput( polydata);
    }
    else{
      normals->SetInput(reader->GetOutput());
    }
  
    // send through normals filter to ensure that triangle orientations
    // are consistent 
    normals->ConsistencyOn();
    normals->SplittingOff();
    normals->AutoOrientNormalsOn();
    
    mesh = normals->GetOutput();
    mesh->Update();
    std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	      << std::endl;

    //Following few lines of code are meant to obtain number of edges
    //from the mesh
    vtkSmartPointer<vtkExtractEdges> extractedEdges = 
      vtkSmartPointer<vtkExtractEdges>::New();
    extractedEdges->SetInput(mesh);
    extractedEdges->Update();
  
    //Number of Cells in vtkExtractEdges = number of edges : Amit
    std::cout << "Number of edges in the mesh= " 
	      << extractedEdges->GetOutput()->GetNumberOfCells()
	      <<std::endl;
  }   


  // create vector of nodes
  double Rcapsid = 1.0;
  
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  
  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;

    if(continueFromNum > 0){
      DeformationNode<3>::Point d;
      DeformationNode<3>::Point temp;
      displacements->GetTuple(a,&(d[0]));
      mesh->GetPoint(a, &(temp[0]));
      x = temp + d;
    }
    else{
      mesh->GetPoint(a, &(x[0]));
      Ravg += tvmet::norm2(x);
    }
    
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);    
    nodes.push_back( n );
    defNodes.push_back( n );
  }

  assert(nodes.size()!=0);

  if(continueFromNum == 0){
    Ravg /= nodes.size();
  }
  cout << "Number of nodes: " <<nodes.size() << endl
       << "Ravg = " << Ravg << endl;

  // read in triangle connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=mesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  if(verbose) cout << "Number of triangles: " <<ntri << endl;
  for (int i = 0; i<ntri; i++){
    assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
    for(int a=0; a<3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
    connectivities.push_back(c);
  }

  std::vector<double> lengthStat;  
  double EdgeLength;
  double stdDevEdgeLen;
  double C0 = 0.0;

  if(continueFromNum == 0){
    // Calculate side lengths average and std dev of the 
    //equilateral triangles
    lengthStat = calcEdgeLenAndStdDev(defNodes,connectivities);  
    EdgeLength = lengthStat[0];
    stdDevEdgeLen = lengthStat[1];
    std::cout<<"Before any relaxation :" << endl
	     <<"   Average triangle edge length = "<< std::setprecision(10)
	     << EdgeLength << endl
	     <<"   Standard deviation = " << std::setprecision(10)
	     << stdDevEdgeLen << endl;
    std::cout.precision(6);
  
    //Rescale the capsid such that triangle edge-lengths are unity
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      x *= 1.0/EdgeLength;
      defNodes[i]->setPoint(x);
      defNodes[i]->setPosition(x);
    }
    
    //Recalculate edge lengths and capsid radius
    lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);  
    EdgeLength = lengthStat[0];

    Ravg = 0.0;
    for(int i=0; i < defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      double tempRadius = tvmet::norm2(x);
      Ravg += tempRadius;
    }
    Ravg /= defNodes.size();

    std::cout<<"Radius of capsid after rescaling = "<< Ravg << endl;
  }

  //Material properties
  if(continueFromNum == 0){
    Rshift = EdgeLength;
  }
  double sigma  = (100/(Rshift*percentStrain))*log(2.0);
  double PotentialSearchRF=1.2*Rshift; 
  double springConstant = 2*sigma*sigma*epsilon;
  double pressure = 12*sigma*epsilon
    *(exp(-2*sigma*Rshift)- exp(-sigma*Rshift)
      + exp(-1.46410*sigma*Rshift) - exp(-0.7321*sigma*Rshift))
    /(3*Ravg*Ravg);

  if(pressure < 0.0){
    pressure = pressure*(-1);
  }

  double fracturePressure = (3.82)*sigma*epsilon/(Rshift*Rshift);
  std::cout<<"Fracture Pressure = "<< fracturePressure << endl
	   <<"Minimum Pressure = "<< pressure << endl;
  pressure *= pressureFactor;
  std::cout<<"Pressure in use = "<< pressure << endl;

  double gamma = gamma_inp;
  double Y = 2.0/sqrt(3)*springConstant; // Young's modulus
  double nu = 1.0/3.0;
  double KC = Y*Ravg*Ravg/gamma; 
  double KG = -2*(1-nu)*KC; // Gaussian modulus
  C0 = 0.0;
  double ARtol = 1.5;
  int quadOrder = 2;

  int m=5;
  int maxIter=1e5;
  double factr=1.0e+1;
  double pgtol=1.0e-7;
  int iprint = 1;

  std::stringstream sstm;
  string fname = modelName;
  string rName;
  string actualFile;

  typedef FVK MaterialType;	
  typedef LoopShellBody<MaterialType> LSB;
  typedef LoopShell<MaterialType> LS;
  MaterialType bending(KC,KG,C0,0.0,0.0);

  if(harmonicRelaxNeeded){

    //****** Relax the initial mesh using harmonic potential ****** //

    LSB * bd1 = new LSB(bending, connectivities, nodes, quadOrder, 
			pressure, 0.0,0.0,1.0e4,1.0e6,1.0e4,
			multiplier,noConstraint,noConstraint);
    bd1->setOutput(paraview);
    SpringPotential SpringMat(springConstant, Rshift);
    PotentialBody * SpringBody = new 
      PotentialBody(&SpringMat, defNodes, PotentialSearchRF);

    //Create Model
    Model::BodyContainer bdc1;
    bdc1.push_back(SpringBody);
    bdc1.push_back(bd1);    
    Model model1(bdc1,nodes);
     
    std::cout<< "Spring constant: " << springConstant << endl;
    std::cout<< "Relaxing the mesh using harmonic potential..."<< endl;

    Lbfgsb solver(model1.dof(), m, factr, pgtol, iprint, 1e5 );
    solver.solve( &model1 );
    
    std::cout<<"Harmonic potential relaxation complete." << endl;

    //Print to VTK file
    sstm << fname <<"-relaxed-harmonic";
    rName = sstm.str();
    model1.print(rName);
    sstm <<"-bd1.vtk";
    actualFile = sstm.str();
    sstm.str("");
    sstm.clear();
    sstm << fname <<"-relaxed-harmonic.vtk";
    rName = sstm.str();
    std::rename(actualFile.c_str(),rName.c_str());      
    sstm.str("");
    sstm.clear();

    //Release allocated memory
    delete bd1;
    delete SpringBody;

    //Recalculate edge lengths and dependent quantities
    lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);  
    EdgeLength = lengthStat[0];
    Rshift = EdgeLength;
    sigma  = (100/(Rshift*percentStrain))*log(2.0);
    springConstant = 2*sigma*sigma*epsilon;
    std::cout<<"After relaxing with harmonic potential: "<< endl
	     <<"   Average triangle edge length = "<< std::setprecision(10)
	     << EdgeLength << endl
	     <<"   Standard deviation = " << std::setprecision(10)
	     << stdDevEdgeLen << endl;
    std::cout.precision(6); 
  }

  //********************************** Actual Indentation **********************//

  LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
  bd->setOutput(paraview);

  // If we want pre-stress removed, then reset the reference
  // configuration for stretching, and reset the spontaneous curvature
  // for bending.
  //
  if( prestressFlag == "no" /*|| prestressFlag == "spherical"*/ ) {
    bd->resetReference();
    for(Body::ElementIterator e=bd->elements().begin(); e!=bd->elements().end(); e++) {
      LoopShell<MaterialType> * lse = (LoopShell<MaterialType>*)(*e);
      for(LoopShell<MaterialType>::QuadPointIterator p=lse->quadraturePoints().begin();
	  p!=lse->quadraturePoints().end(); p++) {
	p->material.setSpontaneousCurvature(2.0*( p->material.meanCurvature() ) );
      }
    }
    bd->SetRefConfiguration(EdgeLength);
  }

  std::cout<< "Morse  potential parameters:" << endl
	   << "sigma = " << sigma << endl	   
	   << "PotentialSearchRF = " << PotentialSearchRF << endl
	   << "epsilon = " << epsilon << endl
	   << "Is Remesh On? " << remesh << endl
	   << "ARtol = " << ARtol << endl;
  
  // Protein body implemented using Morse potential body
  Morse Mat(epsilon, sigma, Rshift);

  // Then initialize potential body
  PotentialBody * PrBody = new PotentialBody(&Mat,defNodes,PotentialSearchRF);  
  PrBody->compute(true, false, false);
  std::cout << "Initial protein body energy = " << PrBody->energy() << endl;    

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(PrBody);
  bdc.push_back(bd);
  
  Model model(bdc,nodes);

  Lbfgsb solver2(model.dof(), m, factr, pgtol, iprint, 1e5);

  if(continueFromNum == 0){

    std::cout << "Relaxing shape for gamma = " << gamma_inp << std::endl;
  
    for(int n=0; n<nodes.size(); n++) {
      for(int i=0; i<nodes[n]->dof(); i++) nodes[n]->setForce(i,0.0);
    }
    
    for(int b=0; b<bdc.size(); b++) {
      std::cout << "bdc[" << b << "]->compute()" << std::endl;
      bdc[b]->compute(true,true,false);    
    }
    
    std::cout << "Initial Shape." << std::endl
	      << "Energy = " << solver2.function() << std::endl;
    
    fname = modelName;
    fname += ".initial";
    model.print(fname);
    actualFile = fname + "-bd1.vtk";
    std::rename(actualFile.c_str(),fname.c_str());
    
    // relax initial shape;
    solver2.solve(&model);
    
    std::cout << "Shape relaxed." << std::endl
	      << "Energy = " << solver2.function() << std::endl;
    
    fname = modelName;
    fname += ".relaxed";
    model.print(fname);
    actualFile = fname + "-bd1.vtk";
    std::rename(actualFile.c_str(),fname.c_str());

    //Calculate centre of sphere as average of position vectors of all nodes.
    tvmet::Vector<double,3> Xavg(0.0);
    for ( int i = 0; i<defNodes.size(); i++){
      Xavg += defNodes[i]->point();
    }
    Xavg /= defNodes.size();
    
    //We will calculate radius using the quadrature points
    LSB::FeElementContainer elements = bd->shells();
    std::vector<double> qpRadius(elements.size(),0.0);
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int e=0; e < elements.size();e++){
      
      const LS::NodeContainer eleNodes = elements[e]->nodes();      
      LS::QuadPointContainer quadPoints = elements[e]->quadraturePoints();
      
      for(LS::ConstQuadPointIterator quadPoint = quadPoints.begin();
	  quadPoint != quadPoints.end(); ++quadPoint){
	
	LoopShellShape s = (*quadPoint).shape;
	const LoopShellShape::FunctionArray fn = s.functions();
	tvmet::Vector<double,3> Xq(0.0);
	
	for(int i=0;i<fn.size();i++){
	  Xq += tvmet::mul(eleNodes[i]->point(),fn(i));
	}
	
	double qpR = tvmet::norm2(Xq-Xavg);
	qpRadius[e] = qpR;
      }
    }
    
    Ravg = 0.0;
    for ( int i = 0; i < qpRadius.size(); i++){
      Ravg += qpRadius[i];
    }
    Ravg /= qpRadius.size();
    std::cout<<"Radius of capsid after relaxation = "<< Ravg << endl;
  
    double dRavg2 = 0.0;
    for ( int i = 0; i<qpRadius.size(); i++){
      double dR =  qpRadius[i]-Ravg; 
      dRavg2 += dR*dR;
    }
    dRavg2 /= qpRadius.size();
    
    double asphericity = dRavg2/(Ravg*Ravg);
    double gammaCalc = Y*Ravg*Ravg/KC;

    std::cout << "Effective 2D Young's modulus = " << Y << endl
	      << "Effective FVK number = " << gammaCalc << endl
	      << "Asphericity = " << asphericity << endl;
    
    // find top and bottom of capsid
    Zmin=std::numeric_limits<double>::max();
    Zmax=-std::numeric_limits<double>::max();
    
    for (int a=0; a<defNodes.size(); a++){
      double Z = defNodes[a]->getPoint(2);
      Zmin = std::min(Zmin, Z);
      Zmax = std::max(Zmax, Z);
    }
    
    std::cout<< "Zmax = " << Zmax << std::endl
	     << "Zmin = " << Zmin;

    xc = 0.0, 0.0, Zmax+afmR;
  }
  
  std::cout << "AFM Indenter radius =" << afmR << std::endl
	    << "AFM Indenter center = (" << xc[0] << "," 
	    << xc[1] << "," << xc[2] << ")" << std::endl;

  std::cout << "Compressing capsid." << std::endl;


  double friction = friction_inp;

  double k_AL = 1.0e2;
  RigidHemisphereAL * afm 
    = new RigidHemisphereAL( defNodes, k_AL, afmR, xc, friction );
  afm->updateContact();
  
  bd->pushBack( afm );
  std::cout << "Added afm to body." << std::endl;

  bool up=true;
  
  if(continueFromNum ==0 ) Z_glass = Zmin;

  RigidPlateAL* glass  = new RigidPlateAL(defNodes, k_AL, Z_glass, up, friction);
  bd->pushBack( glass );

  const double originalHeight=Zmax-Zmin;
  double Zbegin=Z_glass;
  double Zend=Zmin+1.0*Ravg;
  if(indent_inp > 0.0) {
    Zend = Zbegin+indent_inp*Ravg;
  }
  
  if(continueFromNum == 0){
    dZ = (Zend-Zbegin)/100;  
    if(step_inp > 0.0) {
      dZ = step_inp*Ravg;
    }
  }

  // add some viscosity for regularization
  double minViscosity = 1.0e-6*viscosity_inp;
  double maxViscosity = 1.0e+6*viscosity_inp;
  ViscousRegularizer vr(bd->nodes(), viscosity_inp);
  bd->pushBack( &vr ); 

  // set viscosity parameters
  double targetVelocity = std::abs(dZ);
  double vrTol = 1.0e-10;

  double vrEnergy = vr.energy();
  double bdEnergy = bd->energy();

  string fzName = modelName + ".fz";
  ofstream FvsZ(fzName.c_str());

  double F_prev = -1.0;
  double Z_drop = Zbegin;

  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());

  typedef std::vector<DeformationNode<3>* >::const_iterator ConstNodeIterator; 

  if(continueFromNum == 0){
    model.getField(solver2);
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver2.field(i);
    u_prev = 0.0;
  }
  else{
    for(ConstNodeIterator n=defNodes.begin(); n!=defNodes.end(); n++) {
      double X_ref[3];
      double disp[3];
      mesh_prev->GetPoint((*n)->id(), X_ref);
      displacements_prev->GetTuple((*n)->id(), disp);

      const NodeBase::DofIndexMap & idx = (*n)->index();
      for(int ni=0; ni<(*n)->dof(); ni++){
	x_prev(idx[ni]) = X_ref[ni] + disp[ni];
	u_prev(idx[ni]) = disp[ni];
      }
      
    }
  }
  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%

  int step=continueFromNum;

  // Following variables is for output from LoopShellBody::Remesh()
  uint elementsChanged = 0;

  for(double Z = Zbegin; ;Z+=dZ, step++) {
    
    // initial guess
    if(step==0) {
      // shift capsid up by dZ/2 as an initial guess
      for(int a=0; a<defNodes.size(); a++ ) {
	defNodes[a]->addPoint(2,0.5*dZ);
      }
      model.getField(solver2);
    } else if( std::abs(dZ) > 0 ){ // 
      // add scaled version of previous displacement as an initial
      // guess, but only in loading direction
      for(int i=0; i<model.dof(); i++ ) {
	solver2.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
      }      
      model.putField(solver2);
    }

    // move glass up by dZ
    glass->setZ(Z);


    std::cout << std::endl
	      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	      << std::endl
	      << " Z = " << Z 
	      << " zeta = " << Z-Zbegin 
	      << std::endl;

    //viscosity = std::max(viscosity,viscosity_inp);
    //vr.setViscosity(viscosity_inp);
    int viterMax = 20;
    for(int viter=0; viter<viterMax; viter++) {

      if(viter==viterMax-1) vr.setViscosity(minViscosity);

      blitz::Array<double,1> vSave(model.dof());
      model.getField(solver2);
      for(int i=0; i<model.dof(); i++ ) vSave(i) = solver2.field(i);
      
      if(verbose) std::cout << std::endl 
			    << "VISCOUS ITERATION: " << viter 
			    << "\t viscosity = " << vr.viscosity()
			    << std::endl
			    << std::endl;

      // update contact
      glass->updateContact();
      afm->updateContact();


      model.computeAndAssemble(solver2, false, true, false);
      model.print("contact");

      
      solver2.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bd->energy();

      if(verbose) {
	std::cout << "ENERGY:" << std::endl
		  << "viscous energy = " << vrEnergy << std::endl
		  << "   body energy = " << bdEnergy << std::endl
		  << "  total energy = " << solver2.function() << std::endl
		  << std::endl;
      }

      // update contact
      if(verbose) {
	std::cout << "CONTACT:" << std::endl;
	std::cout << "       top active  = " << afm->active() << std::endl
		  << "    bottom active  = " << glass->active() << std::endl
		  << "        top force  = " << afm->FZ() << std::endl 
		  << "     bottom force  = " << glass->FZ() << std::endl
		  << "  top penetration  = " << afm->penetration() << std::endl
		  << "bottom penetration = " << glass->penetration() << std::endl;
      }
	
      if(verbose) {
	std::cout << "VISCOSITY: " << std::endl
		  << "          velocity = " << vr.velocity() << std::endl
		  << "   target velocity = " << targetVelocity << std::endl
		  << " updated viscosity = " << vr.viscosity() << std::endl
		  << std::endl;
      }

      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < std::abs(vrTol*bdEnergy) && solver2.projectedGradientNorm()<=pgtol &&
	 std::max( afm->penetration(), glass->penetration() ) < 1.0e-2*std::abs(dZ) ) {
	// viscous energy is small enough; exit
	break;
      }
    }

    double height = Zmax-Z;
    // add up forces on top and bottom
    double Ztop = Zmax;
    double Zbot = Zmax - height;       

    if(verbose) {
      std::cout << "Contact and Viscous energy converged." << std::endl
		<< std::endl
		<< "        height = " << height << std::endl
		<< "   indentation = " << originalHeight-height << std::endl
		<< std::endl;
    }

    //
    // Do we want to keep this solution or back up and try again?
    //

    //
    // Keeping solution
    //

    // Compute displacement
    for(int i=0; i<model.dof(); i++ ) {
      u_prev(i) = ( solver2.field(i) - x_prev(i) )/std::abs(dZ);
    }

    // If displacement was bigger than dZ, some buckling event must
    // have occurred and it's better not to make a continuation
    // attempt.
    if( max( abs( u_prev ) ) > 1.0 ) u_prev = 0.0;

    // Save current (successful) state as previous
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver2.field(i);

    if ( unload && Z >= Zend ) { // reached max indentation, now unload by reversing dZ
      dZ = -dZ;
    } 
    F_prev = afm->FZ();

    //This is where we remesh before next indentation step

    // Do this only for the bending body since LoopShellBody is the
    // only one that has remeshing implemented

    if(remesh) {     
      elementsChanged = bd->Remesh(ARtol,bending,quadOrder);

      //Print out the number of elements that changed due to remeshing
      if(elementsChanged > 0){
	std::cout<<"Number of elements that changed after remeshing = "
		 <<elementsChanged<<"."<<std::endl;

	//If some elements have changed then we need to reset the
	//reference configuration with average side lengths
	bd->SetRefConfiguration(EdgeLength);

	//We also need to recompute the neighbors for PotentialBody
	PrBody->recomputeNeighbors(PotentialSearchRF);
	
	//Relax again after remeshing
	solver2.solve( &model );

      }

    }// Remeshing ends here

    //*********** BEGIN PRINTING OUTPUT (and log) FILES ***********//

    FvsZ << std::setw( 24 ) << std::setprecision(16) 
	 << originalHeight-height
	 << std::setw( 24 ) << std::setprecision(16) 
	 << glass->FZ()
	 << std::setw( 24 ) << std::setprecision(16) 
	 << afm->FZ()
	 << std::setw( 10 ) 
	 << step
	 << std::endl;
   
    for(int b=0; b<bdc.size(); b++) {
      char name[100];      
      sprintf(name,"%s-body%d-step%04d",modelName.c_str(),b,step);
      bdc[b]->printParaview(name);
      //We will append Caspsomer POINT_DATA to the vtk output file
      //printed by printParaview(), if such a file exists
      sprintf(name,"%s-body%d-step%04d.vtk",modelName.c_str(),b,step);
      if (ifstream(name)){      
	insertValenceInVtk(name,mesh);
      }
    }
    //************* END PRINTING OUTPUT FILES **************//

    // check if we are done
    if( unload && Z+dZ < Zbegin ) 
      break;
    else if( !unload && Z+dZ > Zend ) 
      break;
    
  }// Indentation Loop Ends

  FvsZ.close();

  std::cout << "Indentation complete." << std::endl; 
 
  t2=clock();
  float diff ((float)t2-(float)t1);
  std::cout<<"Total execution time: "<<diff/CLOCKS_PER_SEC
	   <<" seconds"<<std::endl;
  return 0;

}

///////////////////////// END OF MAIN FUNCTION //////////////////////////////////
/////////////////////////                     //////////////////////////////////

///////////////////////// INSERTVALENCEINVTK BEGINS ///////////////////////////
/////////////////////////                           //////////////////////////

//The method insertValenceInVtk() inserts valence information in a vtk file
//The calling method must ensure that 'fileName' exists and is a valid vtk file.
//MUST use the extension '.vtk' in 'fileName'.
//'mesh' is a pointer to a vtkPolyData

void insertValenceInVtk(const std::string fileName, vtkSmartPointer<vtkPolyData> mesh){
  ofstream * appendTo;

  //Check that the file exists
  assert(ifstream(fileName.c_str()));

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName(fileName.c_str());
  
  vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();  
  ds->Update();
  
  //The following vtkUnsignedIntArray will be used to store the
  //number of CELLS in the mesh that share the POINT denoted by the
  //index of the vector
  vtkSmartPointer<vtkUnsignedIntArray> countPointCells = 
    vtkSmartPointer<vtkUnsignedIntArray>::New();
  countPointCells->SetNumberOfValues(mesh->GetNumberOfPoints());
  countPointCells->SetName("Valence");
  
  //cellIds will be used to temporarily hold the CELLS that use a
  //point specified by a point id.
  vtkSmartPointer<vtkIdList> cellIds = 
    vtkSmartPointer<vtkIdList>::New();
  
  //We will use a vtkPolyDataWriter to write our modified output
  //files that will have Capsomer information as well
  vtkSmartPointer<vtkDataWriter> writer = 
    vtkSmartPointer<vtkDataWriter>::New();
  
  if(ds->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){
    vtkSmartPointer<vtkUnstructuredGrid> usg 
      = vtkUnstructuredGrid::SafeDownCast(ds);
    usg->BuildLinks();
    for(int p=0; p<ds->GetNumberOfPoints(); p++){
      usg->GetPointCells(p,cellIds);
      countPointCells->SetValue(p,cellIds->GetNumberOfIds());
      cellIds->Reset();
    }      
  }                  
  else if(ds->GetDataObjectType() == VTK_POLY_DATA){
    vtkSmartPointer<vtkPolyData> pd = vtkPolyData::SafeDownCast(ds);
    pd->BuildLinks();
    for(int p=0; p<ds->GetNumberOfPoints(); p++){
      pd->GetPointCells(p,cellIds);
      countPointCells->SetValue(p,cellIds->GetNumberOfIds());
      cellIds->Reset();
    }	
  }
  ds->GetFieldData()->AddArray(countPointCells);
  appendTo = new ofstream(fileName.c_str(),ofstream::app);
  writer->WriteFieldData(appendTo,ds->GetFieldData());
  appendTo->close();
}

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                    CALCEDGELENANDSTDDEV BEGINS                        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
/*
  Calculates average edge lengths of triangles in the mesh and the
  standard deviation in the edge lengths.
*/

std::vector<double> calcEdgeLenAndStdDev
(std::vector< DeformationNode<3>* > defNodes, 
 vector< tvmet::Vector<int,3> > connectivities){

  double EdgeLength = 0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> cm(3);
    for(int j=0; j<3; j++) cm[j]=connectivities[i](j);
    // Edge vectors in current config.
    tvmet::Vector<double,3> 
      e31(defNodes[cm[0]]->point()-defNodes[cm[2]]->point()), 
      e32(defNodes[cm[1]]->point()-defNodes[cm[2]]->point()),
      e12(defNodes[cm[1]]->point()-defNodes[cm[0]]->point()),
      eCent(defNodes[cm[2]]->point());
    // Compute average edge length for each triangle
    double temp = 
      (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;

#pragma omp atomic
    EdgeLength += temp;
  }
  EdgeLength /= connectivities.size();  

  // Calculate the standard deviation of side lengths of the
  // equilateral triangles
  double stdDevEdgeLen = 0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i=0; i<connectivities.size(); i++) {
    std::vector<int> cm(3);
    for(int j=0; j<3; j++) cm[j]=connectivities[i](j);
    // Edge vectors in current config.
    tvmet::Vector<double,3> 
      e31(defNodes[cm[0]]->point()-defNodes[cm[2]]->point()), 
      e32(defNodes[cm[1]]->point()-defNodes[cm[2]]->point()),
      e12(defNodes[cm[1]]->point()-defNodes[cm[0]]->point()),
      eCent(defNodes[cm[2]]->point());
    double temp = std::pow(tvmet::norm2(e31) - EdgeLength,2.0) +
      std::pow(tvmet::norm2(e32) - EdgeLength,2.0) +
      std::pow(tvmet::norm2(e12) - EdgeLength,2.0);
    
#pragma omp atomic
    stdDevEdgeLen += temp;
  }

  stdDevEdgeLen /= connectivities.size();
  stdDevEdgeLen = sqrt(stdDevEdgeLen);

  std::vector<double> result;
  result.push_back(EdgeLength);
  result.push_back(stdDevEdgeLen);
  return result;
}

