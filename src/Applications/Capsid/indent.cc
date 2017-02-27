#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>

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

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void insertValenceInVtk(const std::string fileName,\
			vtkSmartPointer<vtkPolyData> mesh);

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
  // double remesh_inp=0.0; //
  bool remesh = false;
  string prestressFlag = "yes";
  bool CST=false;
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
      //     case 'n':
      //       maxIter_inp = std::atof(optarg);
      //       std::cout << "max iterations: " << maxIter_inp << std::endl;
      //       break;
    case 'p' :
      //       if( std::string(optarg) == std::string("no") ) 
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
    case 'e':
      if( std::string(optarg) == std::string("CST") ){
	CST = true;
	std::cout << "Using mixed formulation with CST elements for stretching."
		  << std::endl;
      }
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

  string inputFileName = modelName + ".vtk";
  vtkDataSetReader * reader = vtkDataSetReader::New();
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
  vtkSmartPointer<vtkPolyData> mesh = normals->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  //Following few lines of code are meant to obtain number of edges
  //from the mesh
  vtkSmartPointer<vtkExtractEdges> extractedEdges = 
    vtkSmartPointer<vtkExtractEdges>::New();
  extractedEdges->SetInput(mesh);
  extractedEdges->Update();
  //Uncommenting the next line fetches the line entities, if you want 
  //them : Amit
  //vtkCellArray * lines = extractedEdges()->GetOutput()->GetLines();

  //Number of Cells in vtkExtractEdges = number of edges : Amit
  std::cout << "Number of edges in the mesh= " 
	    << extractedEdges->GetOutput()->GetNumberOfCells()
	    <<std::endl; 


  // create vector of nodes
  double Rcapsid = 1.0;
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  // vector<ProteinNode *> Proteins;
  double Ravg = 0;

  // read in points
  for(int a=0; a<mesh->GetNumberOfPoints(); a++) {
    int id=a;
    DeformationNode<3>::Point x;
    mesh->GetPoint(a, &(x[0]));
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    DeformationNode<3>* n = new DeformationNode<3>(id,idx,x);    
    nodes.push_back( n );
    defNodes.push_back( n );
    //We will add a Protein at every node - Amit
    //ProteinNode * PrNode = new ProteinNode(n);
    //Proteins.push_back(PrNode);
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
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

  double C0 = 0.0;

  // rescale size 
  if( prestressFlag == "spherical" ) {
    // make capsid spherical 

    C0 = - 2.0/Rcapsid; // WSK: use minus sign here consistent with
			// outward pointing surface normals.

    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      double R = norm2(x);
      x *= Rcapsid/R;
      defNodes[i]->setPoint(x);
      defNodes[i]->setPosition(x);
    }
  }
  else {
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      x *= Rcapsid/Ravg;
      defNodes[i]->setPoint(x);
      defNodes[i]->setPosition(x);
    }
  }

  double Y = sqrt(gamma_inp);
  double KC = 1.0/Y;
  

  // Introduce a numerical scaling factor to make energy and forces
  // large enough to avoid issues with machine precision.  Later
  // divide forces by factor when printing to file.
  double scalingFactor = 1.0;//e6;
  Y *= scalingFactor;
  KC *= scalingFactor;

  //Amit: Set nu=0 and Yb=0 if you want LoopShellBody not to handle
  //stretching anymore.
  double nu = 1.0/3.0;
  double Yb=Y;
  double KG = -2.0*(1.0-nu)*KC;

  if(verbose) 
    std::cout << " Y: " << Y << std::endl
	      << " nu: " << nu << std::endl
	      << " KC: " << KC << std::endl
	      << " KG: " << KG << std::endl
	      << " C0: " << C0 << std::endl;

  // create Body
  int quadOrder = 2;

  typedef FVK MaterialType;

  MaterialType bending( KC, KG, C0, Yb, nu );
	
  typedef LoopShellBody<MaterialType> LSB;
  LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
  
  bd->setOutput(paraview);

  double EquilateralEdgeLength = 0.0;

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
    EquilateralEdgeLength += 
      (tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;
  }
  EquilateralEdgeLength /= connectivities.size();
  std::cout<<"Average equilateral triangle edge length:"<<EquilateralEdgeLength<<endl;
  
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
    bd->SetRefConfiguration(EquilateralEdgeLength);
  }

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(bd);

  double ARtol = 1.1;    
  
  Model model(bdc,nodes);

  int m=5;
  int maxIter=500;//1000;//model.dof();
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = 0;
  int maxIter_inp=0;
  ifstream lbfgsbinp("lbfgsb.inp");
  
  lbfgsbinp >> iprint >> factr >> pgtol >> m >> maxIter_inp;
  if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input maxIter: " << maxIter_inp << std::endl;
  maxIter = std::max(maxIter,maxIter_inp);

  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );//(true);

  std::cout << "Relaxing shape for gamma = " << gamma_inp << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  for(int n=0; n<nodes.size(); n++) {
    for(int i=0; i<nodes[n]->dof(); i++) nodes[n]->setForce(i,0.0);
  }
    
  for(int b=0; b<bdc.size(); b++) {
    std::cout << "bdc[" << b << "]->compute()" << std::endl;
    bdc[b]->compute(true,true,false);    
  }

  std::cout << "Initial Shape." << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  string fname = modelName;
  fname += ".initial";
  model.print(fname);

  // relax initial shape
  solver.solve(&model);

  std::cout << "Shape relaxed." << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  fname = modelName;
  fname += ".relaxed";
  model.print(fname);

  tvmet::Vector<double,3> Xavg(0.0);
  for ( int i = 0; i<defNodes.size(); i++){
    Xavg += defNodes[i]->point();
  }
  Xavg /= nodes.size();
  Ravg = 0.0;
  for ( int i = 0; i<defNodes.size(); i++){
    Ravg += tvmet::norm2( defNodes[i]->point() - Xavg );
  }
  Ravg /= nodes.size();
  
  double dRavg2 = 0.0;
  for ( int i = 0; i<defNodes.size(); i++){
    double dR = (tvmet::norm2( defNodes[i]->point() - Xavg) - Ravg); 
    dRavg2 += dR*dR;
  }
  dRavg2 /= nodes.size();
  
  double gammaCalc = Y*Ravg*Ravg/KC;
  double asphericity = dRavg2/(Ravg*Ravg);
  
  std::cout << "gamma = " << gammaCalc << endl
	    << "asphericity = " << asphericity << endl;

  std::cout << "Compressing capsid." << std::endl;

  // find top and bottom of capsid
  double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  double Zavg = 0.0;
  for (int a=0; a<defNodes.size(); a++){
    double Z = defNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
    Zmax = std::max(Zmax, Z);
    Zavg += Z;
  }
  Zavg /= defNodes.size();

  // create indentor and plate
  double afmR = 1.0;
  tvmet::Vector<double,3> xc; xc = 0.0, 0.0, Zmax+afmR; 

  double friction = friction_inp;

  double k_AL = 1.0e2*scalingFactor;
  RigidHemisphereAL * afm 
    = new RigidHemisphereAL( defNodes, k_AL, afmR, xc, friction );
  afm->updateContact();
  
  bd->pushBack( afm );
  std::cout << "Added afm to body." << std::endl;

  bool up=true;
  RigidPlateAL* glass  = new RigidPlateAL(defNodes, k_AL, Zmin, up, friction);
  bd->pushBack( glass );

  const double originalHeight=Zmax-Zmin;
  double Zbegin=Zmin;
  double Zend=Zmin+1.0*Ravg;
  double dZ = (Zend-Zbegin)/100;
  if(indent_inp > 0.0) {
    Zend = Zbegin+indent_inp;
  }
  if(step_inp > 0.0) {
    dZ = step_inp;
  }

  // add some viscosity for regularization
  double minViscosity = 1.0e-6*viscosity_inp;
  double maxViscosity = 1.0e+6*viscosity_inp;
  ViscousRegularizer vr(bd->nodes(), viscosity_inp);
  bd->pushBack( &vr ); 

  // set viscosity parameters
  double targetVelocity = dZ;
  double vrTol = 1.0e-10;

  double vrEnergy = vr.energy();
  double bdEnergy = bd->energy();

  string fzName = modelName + ".fz";
  ofstream FvsZ(fzName.c_str());

  double F_prev = -1.0;
  double Z_drop = Zbegin;
  blitz::Array<double,1> x_prev(model.dof());
  blitz::Array<double,1> u_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);
  u_prev = 0.0;

  // %%%%%%%%%%%%%%%%%%%%%%
  // Begin indentation loop
  // %%%%%%%%%%%%%%%%%%%%%%

  int step=0;

  // Following variables is for output from LoopShellBody::Remesh()
  uint elementsChanged = 0;

  for(double Z = Zbegin; /*Z<Zend+0.5*dZ*//*Z>=Zbegin*/; Z+=dZ, step++) {
    
    // initial guess
    if(step==1) {
      // shift capsid up by dZ/2 as an initial guess
      for(int a=0; a<defNodes.size(); a++ ) {
	defNodes[a]->addPoint(2,0.5*dZ);
      }
      model.getField(solver);
    } else if( dZ > 0 ){ // 
      // add scaled version of previous displacement as an initial
      // guess, but only in loading direction
      for(int i=0; i<model.dof(); i++ ) {
	solver.field(i) = x_prev(i) + 0.99*dZ*u_prev(i);
      }      
      model.putField(solver);
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
      model.getField(solver);
      for(int i=0; i<model.dof(); i++ ) vSave(i) = solver.field(i);
      
      if(verbose) std::cout << std::endl 
			    << "VISCOUS ITERATION: " << viter 
			    << "\t viscosity = " << vr.viscosity()
			    << std::endl
			    << std::endl;

      // update contact
      glass->updateContact();
      afm->updateContact();


      model.computeAndAssemble(solver, false, true, false);
      model.print("contact");

      
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bd->energy();

      if(verbose) {
	std::cout << "ENERGY:" << std::endl
		  << "viscous energy = " << vrEnergy << std::endl
		  << "   body energy = " << bdEnergy << std::endl
		  << "  total energy = " << solver.function() << std::endl
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

      if(vrEnergy < std::abs(vrTol*bdEnergy) && solver.projectedGradientNorm()<=pgtol &&
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
      u_prev(i) = ( solver.field(i) - x_prev(i) )/std::abs(dZ);
    }

    // If displacement was bigger than dZ, some buckling event must
    // have occurred and it's better not to make a continuation
    // attempt.
    if( max( abs( u_prev ) ) > 1.0 ) u_prev = 0.0;

    // Save current (successful) state as previous
    for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);

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
	bd->SetRefConfiguration(EquilateralEdgeLength);
	
      }
    }// Remeshing ends here

    //*********** BEGIN PRINTING OUTPUT (and log) FILES ***********//

    FvsZ << std::setw( 24 ) << std::setprecision(16) 
	 << originalHeight-height
	 << std::setw( 24 ) << std::setprecision(16) 
	 << glass->FZ()/scalingFactor
	 << std::setw( 24 ) << std::setprecision(16) 
	 << afm->FZ()/scalingFactor
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



