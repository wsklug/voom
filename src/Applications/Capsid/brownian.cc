#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <iomanip>
#include <limits>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Quadrature.h"
#include "TriangleQuadrature.h"

#include <vtkPolyData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkGeometryFilter.h>
#include <vtkSetGet.h>

#include "Morse.h"
#include "SpringPotential.h"
#include "PotentialBody.h"
#include "BrownianKick.h"
#include "ViscousRegularizer.h"

#include "HelperFunctions.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

typedef DeformationNode<3> Node; // nickname for mechanics nodes
typedef std::vector< Node* > NodeContainer;

int main(int argc, char* argv[])
{
  clock_t t1,t2,t3;
  t1=clock();
  if ( argc != 2 ){
    cout<<"usage: "<< argv[0] <<" <filename>\n";
    return -1;
  }

#ifdef _OPENMP
  std::cout<< "************* PARALLELIZATION USING OPENMP ****************"
	   << endl << endl;
#endif

  string inputFileName = argv[1];

  //For Morse material 
  double epsilon;
  double FVKnum;
  double percentStrain;
  double pressureFactor;

  //Numerical viscosity input parameter
  //double temperature;
  //double viscosity_inp;
  double dt = 9.76e-4; //time step
  int viterMax;

  double Rshift;
  bool rescale;
  int continueFromNum = 1;
  int nameSuffix = 0;
  int firstFileNum;
  bool serialFlag;

  //Read epsilon and percentStrain from input file. percentStrain is
  //calculated so as to set the inflection point of Morse potential
  //at a fixed distance relative to the equilibrium separation
  //e.g. 1.1*R_eq, 1.5*R_eq etc.
  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  string temp;
  miscInpFile >> temp >> epsilon
	      >> temp >> FVKnum
	      >> temp >> Rshift
	      >> temp >> percentStrain
	      >> temp >> pressureFactor
    	      >> temp >> dt
	      >> temp >> viterMax
	      >> temp >> rescale
	      >> temp >> continueFromNum
    	      >> temp >> serialFlag;

  miscInpFile.close();

  vtkSmartPointer<vtkDataSetReader> reader = 
    vtkSmartPointer<vtkDataSetReader>::New();

  std::stringstream sstm;
  
  if(continueFromNum > 1){
    nameSuffix = continueFromNum - 1;
    firstFileNum = nameSuffix;   
    temp = inputFileName.substr(0,inputFileName.find("."));
    sstm << temp <<"-relaxed-"<<firstFileNum-1<<".vtk";
    temp = sstm.str();

    reader->SetFileName(temp.c_str());

    sstm.str("");
    sstm.clear();
    std::cout<<"Continuing from step number "<< continueFromNum 
	     <<"."<< endl <<"Using input file "
	     << temp << endl <<"Rshift = "<< Rshift 
	     << endl << endl;
  }
  else{
    reader->SetFileName( inputFileName.c_str() );
  }

  //We will use this object, shortly, to ensure consistent triangle
  //orientations
  vtkSmartPointer<vtkPolyDataNormals> normals = 
    vtkSmartPointer<vtkPolyDataNormals>::New();

  //We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput() If
  //our input vtk file has vtkUnstructuredGridData instead of
  //vtkPolyData then we need to convert it using vtkGeometryFilter
  vtkSmartPointer<vtkDataSet> ds = reader->GetOutput();
  reader->Update();
  
  vtkSmartPointer<vtkGeometryFilter> geometryFilter = 
      vtkSmartPointer<vtkGeometryFilter>::New();
  
  if(ds->GetDataObjectType() == VTK_UNSTRUCTURED_GRID){      
    geometryFilter->SetInputConnection(reader->GetOutputPort());
    normals->SetInputConnection(geometryFilter->GetOutputPort());
  }
  else{
    normals->SetInputConnection(reader->GetOutputPort());
  }
  
  // send through normals filter to ensure that triangle orientations
  // are consistent 
  normals->ComputeCellNormalsOn();
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  vtkSmartPointer<vtkPolyData> mesh = normals->GetOutput();
  normals->Update();
  
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  // create vector of nodes

  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< DeformationNode<3>* > defNodes;
  double Ravg = 0.0;

  // read in points. DON'T FORGET TO ADD DISPLACEMENTS FROM INPUT FILE IF ANY
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
  }
  assert(nodes.size()!=0);
  Ravg /= nodes.size();
  cout <<"Number of nodes: "<<nodes.size() << endl
       <<"Initial radius: "<<Ravg<<endl;

  // read in triangle connectivities
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=mesh->GetNumberOfCells();
  connectivities.reserve(ntri);
  std::cout << "Number of triangles: " <<ntri << endl;

  for (int i = 0; i<ntri; i++){
    assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
    for(int a=0; a<3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
    connectivities.push_back(c);
  }

  // Calculate side lengths average and std dev of the 
  //equilateral triangles
  std::vector<double> lengthStat = 
    calcEdgeLenAndStdDev(defNodes,connectivities);  
  double EdgeLength = lengthStat[0];
  double stdDevEdgeLen = lengthStat[1];
  std::cout<<"Before any relaxation :" << endl
	   <<"   Average triangle edge length = "<< std::setprecision(10)
	   << EdgeLength << endl
	   <<"   Standard deviation = " << std::setprecision(10)
	   << stdDevEdgeLen << endl;
  std::cout.precision(6);

  if(rescale){
    // Rescale size of the capsid by the average equilateral edge length
    for(int i=0; i<defNodes.size(); i++) {
      DeformationNode<3>::Point x;
      x = defNodes[i]->point();
      x *= 1.0/EdgeLength;
      defNodes[i]->setPoint(x);
      defNodes[i]->setPosition(x);
    }
    //Recalculate edge lengths and dependent quantities
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

  if(continueFromNum == 1){
    Rshift = EdgeLength;
  }

  string fname = inputFileName.substr(0,inputFileName.find("."));
  string iName;
  string rName;
  string actualFile;

  ofstream myfile;
  std::string dataOutputFile = "BrownianRelax.dat";

  //If the file already exists then open it in append mode
  if(!ifstream(dataOutputFile.c_str())){
    myfile.open(dataOutputFile.c_str());
    myfile << setw(5) << "#Step"<<"\t" 
	   << setw(9) << "Ravg" <<"\t"
	   << setw(8) << "asphericity" <<"\t" 
	   << setw(8) << "FVK" <<"\t"
	   << setw(8) << "BendEnergy" <<"\t"
	   << setw(8) << "StretchEnergy" <<"\t" 
	   << setw(8) << "BrownEnergy" <<"\t"
	   << setw(8) << "ViscousEnergy" <<"\t"
	   << setw(8) << "Total Functional" <<"\t"
	   << endl;
    myfile<< showpoint;
  }
  else{
    myfile.open(dataOutputFile.c_str(),ofstream::app);
  }
  
  //Parameters for the l-BFGS solver
  int m=5;

  int maxIter2 = 1e5;
  int maxIter1=1e5;

  double factr=1.0e+1;
  double pgtol=1e-7;
  int iprint = 100;
  Lbfgsb solver1(3*nodes.size(), m, factr, pgtol, iprint, maxIter1 );
  Lbfgsb solver2(3*nodes.size(), m, factr, pgtol, iprint, maxIter2 );

  typedef FVK MaterialType;
  typedef LoopShellBody<MaterialType> LSB;
  typedef LoopShell<MaterialType> LS;

  //****************  Protein body parameters ****************//

  double sigma  = (100/(Rshift*percentStrain))*log(2.0);
  //double PotentialSearchRF = 1.5*Ravg;
  double PotentialSearchRF=2*Rshift; 
  double springConstant = 2*sigma*sigma*epsilon;
  double pressure = 0.0; 

  bool pressureOn = false;
  if(pressureOn){
    pressure = 12*sigma*epsilon
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
  }

  double diffusionCoeff = 4.0*Rshift*Rshift;
  double Cd = 1.0/diffusionCoeff;
  double viscosity = Cd/dt;

  std::cout << "Viscosity Input Parameters:" << std::endl
	    << " Cd = " << Cd << std::endl
	    << "  D = " << diffusionCoeff << std::endl
	    << " dt = " << dt << std::endl;

  double gamma = FVKnum;
  double Y = 2.0/sqrt(3)*springConstant; // Young's modulus
  double nu = 1.0/3.0;
  double KC = Y*Ravg*Ravg/gamma; 
  double KG = -2*(1-nu)*KC; // Gaussian modulus
  double C0 = 0.0;
  int quadOrder = 2;

  double vrEnergy;
  double bkEnergy;
  double bdEnergy;
  double PrEnergy;
  double energy = 0.0;

  //The Bodies
  MaterialType bending(KC,KG,C0,0.0,0.0);
  //LSB * bd = new LSB(bending, connectivities, nodes, quadOrder, pressure,
  //		     0.0,0.0,1.0e4,1.0e6,1.0e4,multiplier,noConstraint,noConstraint);

  LSB * bd;

  bool areaConstraintOn = true;
  if(areaConstraintOn){
    bd = new LSB(bending, connectivities, nodes, quadOrder, pressure,
		       0.0,0.0,1.0e4,1.0e6,1.0e4,multiplier,penalty,noConstraint);
    std::cout << "Prescribed Area = "<< bd->prescribedArea()  << std::endl;
  }
  else{
    std::cout<< "********** CONSTRAINTS NOT BEING USED **********"<<std::endl;
    bd = new LSB(bending, connectivities, nodes, quadOrder);
  }

  bd->setOutput(paraview);

  Morse Mat(epsilon,sigma,Rshift);
  PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);
  ViscousRegularizer vr(bd->nodes(), viscosity);
  bd->pushBack(&vr);
  BrownianKick bk(nodes,Cd,diffusionCoeff,dt);
  bd->pushBack(&bk);

  //Create Model
  Model::BodyContainer bdc;
  bdc.push_back(PrBody);
  bdc.push_back(bd);    
  Model model(bdc,nodes);
     
  std::cout<< "Morse potential parameters:" << endl
	   << "sigma = " << sigma <<" epsilon = " << epsilon
	   << " Rshift = "<< Rshift <<endl;

  std::cout << "Pressure = " << pressure << endl
	    << "Capsid radius = "<< Ravg << endl;

  bool checkConsistency = false;
  if(checkConsistency){
    std::cout<< "Checking consistency......"<<std::endl;
    serialFlag? bk.updateSerialKick(): bk.updateParallelKick();
    bd->checkConsistency(true);
    PrBody->checkConsistency(true);
  }


  //***************************  SOLUTION LOOP ***************************//  
  
  for(int viter = continueFromNum-1; viter < viterMax; viter++) {

    std::cout << std::endl 
	      << "VISCOUS ITERATION: " << viter 
	      << "\t viscosity = " << vr.viscosity()
	      << std::endl
	      << std::endl;
      
    serialFlag? bk.updateSerialKick():bk.updateParallelKick();

    //For debugging we have limited the number of solver iterations to
    //100 so that we can see the intermediate results before the
    //solver diverges
    /*
    for(int z=0; z<interimIter; z++){
	
      solver1.solve( &model );
      //Print the files
      sstm << fname <<"-interim-" << nameSuffix << "-"<< z <<".vtk";
      rName = sstm.str();
      char tempName[] = "lbfgsbconv-bd1.vtk";
      int renameBool = std::rename(tempName,rName.c_str());
      if (renameBool != 0){
	model.print(rName);
	sstm.str("");
	sstm.clear();
	sstm << rName <<"-bd1.vtk";
	actualFile = sstm.str();
	std::rename(actualFile.c_str(),rName.c_str());
      }
      sstm.str("");
      sstm.clear(); // Clear state flags
    }
    */

    //********** Print configuration before relaxation************//
    sstm << fname <<"-initial-" << nameSuffix;
    rName = sstm.str();
    model.print(rName);
    sstm <<"-bd1.vtk";
    actualFile = sstm.str();
    sstm.str("");
    sstm.clear();
    sstm << fname <<"-initial-" << nameSuffix <<".vtk";
    rName = sstm.str();
    std::rename(actualFile.c_str(),rName.c_str());
    sstm.str("");
    sstm.clear();

    solver1.solve( &model );

    vrEnergy = vr.energy();
    bkEnergy = bk.energy();
    bdEnergy = bd->energy() - vrEnergy - bkEnergy;
    PrEnergy = PrBody->energy();
    energy = solver1.function();
    
    std::cout << "ENERGY:" << std::endl
	      << "viscous energy = " << vrEnergy << std::endl
	      << "Brownian energy = " << bkEnergy << std::endl
	      << "protein energy = " << PrEnergy << std::endl
	      << "bending energy = " << bdEnergy << std::endl
	      << "  total energy = " << energy   << std::endl
	      << std::endl;
    std::cout << "VISCOSITY: " << std::endl
	      << "          velocity = " << vr.velocity() << std::endl
	      << " updated viscosity = " << vr.viscosity() << std::endl
	      << std::endl;

    //*******************  REMESHING **************************//
    bool remesh = true;
    double ARtol = 1.2;
    uint elementsChanged = 0;
    
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
	vrEnergy = vr.energy();
	bkEnergy = bk.energy();
	bdEnergy = bd->energy() - vrEnergy - bkEnergy;
	PrEnergy = PrBody->energy();
	energy = solver2.function();
      }
    }
    //*********************************************************//
  
    std::cout << "Shape relaxed." << std::endl
	      << "Energy = " << energy << std::endl;

    //Calculate maximum principal strains in all elements
    //bd->calcMaxPrincipalStrains();

    //********** Print relaxed configuration ************//
    sstm << fname <<"-relaxed-" << nameSuffix;
    rName = sstm.str();
    model.print(rName);
    sstm <<"-bd1.vtk";
    actualFile = sstm.str();
    sstm.str("");
    sstm.clear();
    sstm << fname <<"-relaxed-" << nameSuffix <<".vtk";
    rName = sstm.str();
    std::rename(actualFile.c_str(),rName.c_str());
    sstm.str("");
    sstm.clear();
    //****************************************************//

    //Re-calculate triangle edge length mean and deviation
    lengthStat = calcEdgeLenAndStdDev(defNodes, connectivities);
    EdgeLength = lengthStat[0];
    stdDevEdgeLen = lengthStat[1];
    std::cout<<"After relaxing with Morse potential: "<< endl
	     <<"   Average triangle edge length = "<< std::setprecision(10)
	     << EdgeLength << endl
	     <<"   Standard deviation = " << std::setprecision(10)
	     << stdDevEdgeLen << endl;
    std::cout.precision(6);    

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

    myfile<< nameSuffix++<<"\t\t"<< Ravg <<"\t\t"<< asphericity
	  <<"\t\t"<< gamma <<"\t\t" << bdEnergy <<"\t\t"<< PrEnergy
	  <<"\t\t"<< bkEnergy <<"\t\t" << vrEnergy <<"\t\t"
	  << solver1.function()<< endl;
   
    // step forward in "time", relaxing viscous energy & forces 
    vr.step();
  }
   
  //Release the dynamically allocated memory
  delete bd;
  delete PrBody;
  

  myfile.close();
  t2=clock();
  float diff ((float)t2-(float)t1);
  std::cout<<"Solution loop execution time: "<<diff/CLOCKS_PER_SEC
	   <<" seconds"<<std::endl;

  // Post-processing: Manipulating VTK files
  std::vector<std::string> allVTKFiles;
  int numFiles = viterMax;
  allVTKFiles.reserve(numFiles);
  for(int fileNum=0 ; fileNum < numFiles; fileNum++){    
    sstm << fname <<"-relaxed-" << fileNum <<".vtk";
    std::string tempString = sstm.str();
    allVTKFiles.push_back(tempString);
    sstm.str("");
    sstm.clear();
  }
  
  insertValenceInVtk(allVTKFiles);
  writeEdgeStrainVtk(allVTKFiles,Rshift,percentStrain);
  t3=clock();
  diff = ((float)t3-(float)t2);
  std::cout<<"Post-processing execution time: "<<diff/CLOCKS_PER_SEC
  	   <<" seconds"<<std::endl;
  
}
