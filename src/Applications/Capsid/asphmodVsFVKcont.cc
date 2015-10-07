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
#include "Model.h"
#include "Lbfgsb.h"

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

#include "LennardJones.h"
#include "PotentialBody.h"
#include "Utils/PrintingProtein.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  clock_t t1,t2;
  t1=clock();

  if ( argc != 2 ){
    cout<<"usage: "<< argv[0] <<" <filename>\n";
    return -1;
  }

  bool remesh = true;
  string prestressFlag = "yes";

  string inputFileName = argv[1];
  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );

  //We will use this object, shortly, to ensure consistent triangle
  //orientations
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();

  //We have to pass a vtkPolyData to vtkPolyDataNormals::SetInput() If
  //our input vtk file has vtkUnstructuredGridData instead of
  //vtkPolyData then we need to convert it using vtkGeometryFilter
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
  std::cout << "Number of triangles: " <<ntri << endl;

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
  
  // We want variable number of FVK increments in different ranges. So
  // we will read FVK values from a file instead of generating it by
  // code  
  std::ifstream fvkFile("fvkSteps.dat");
  assert(fvkFile);
  std::vector<vector<double> > gammaVec;
  double currFVK,currPrintFlag;  
  while(fvkFile>>currFVK>>currPrintFlag){
    std::vector<double> currLine;
    currLine.push_back(currFVK);
    currLine.push_back(currPrintFlag);
    gammaVec.push_back(currLine);
  }  
  fvkFile.close();
  
  //Uncomment if you want to calculate FVK values by code
  /*
  double expMin = .0; // log_10(1) = 0 -> FVK = 1
  double expMax = 4.0; // log_10(10000) = 4 -> FVK = 10000
  double FVK_steps = 1000;
  
  // Calculate increment to create 1000 evenly spaced exponents from 0
  // to 4
  double expIncr = (expMax-expMin)/FVK_steps;
  double currExponent; // Used to calculate FVK number
  */
  
  // The Young's modulus and bending modulus
  double Y;
  double KC;
  double nu = 1.0/3.0;
  double KG;
  
  //Order of quadrature
  int quadOrder = 2;

  typedef FVK MaterialType;  	
  typedef LoopShellBody<MaterialType> LSB;
  typedef LoopShell<MaterialType> LS;
  LSB * bd;

  // create Model
  Model::BodyContainer bdc;

  //Parameters for the l-BFGS solver
  int m=5;
  int maxIter=1e5;
  double factr=1.0e+1;
  double pgtol=1.0e-8;
  int iprint = 1;  
 
  std::stringstream sstm;
  string fname = inputFileName.substr(0,inputFileName.find("."));
  string iName;
  string rName;

  ofstream myfile;
  myfile.open ("asphVsFVKcont.dat");
  myfile << "#Ravg\tY\tasphericity\tFVK\tAvgStrain"<< endl;
  myfile << showpoint;

  //Read pressure value from input file
  double pressure;
  std::ifstream miscInpFile("miscInp.dat");
  assert(miscInpFile);
  string temp;
  miscInpFile >> temp >> pressure;
  miscInpFile.close();
  std::cout << "Pressure: " << pressure;

  //Output file index number
  int fileNum = 1;
 
  double asphericity;
  Lbfgsb solver(3*nodes.size(), m, factr, pgtol, iprint, maxIter );

  //***************************  FOR_LOOP ***************************//

  //Loop over all values of FVK number and relax the shape
  for(std::vector<vector<double> >::iterator q=gammaVec.begin();
      q!=gammaVec.end();++q){    
    
    //Uncomment if you want to calculate FVK values by code
    /*
    currExponent = expMin + q*expIncr;
    gamma = pow(10.0,currExponent);
    */
    
    double gamma = (*q)[0];
    double currPrintFlag = (*q)[1];

    // Set the Young's modulus and other material properties based on
    // current FVK number
    Y = sqrt(gamma);
    KC = 1.0/Y;
    KG = -2.0*(1.0-nu)*KC;

    //Remove the bending body from previous iteration and add new
    //one
    if (!bdc.empty()){
      bdc.pop_back();
      delete bd;
    }
    MaterialType bending( KC, KG, C0, Y, nu );
    bd = new LSB(bending, connectivities, nodes, quadOrder,pressure,
		 0.0,0.0,1.0e4,1.0e6,1.0e4,multiplier,noConstraint,noConstraint);    
    bd->setOutput(paraview);
    bdc.push_back(bd);
    Model model(bdc,nodes);
    if(q==gammaVec.begin()){// We need to do this only once
      solver.resize(model.dof());
    }

    std::cout << "Relaxing shape for gamma = " << gamma << std::endl
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
    
    //Uncomment the following region if you want to print initial
    //shapes
    /*
    sstm << fname <<"-initial-" << q;
    iName = sstm.str();

    model.print(iName);

    sstm.str("");
    sstm.clear(); // Clear state flags.
    */

    // relax initial shape
    solver.solve(&model);

    std::cout << "Shape relaxed." << std::endl
	      << "Energy = " << solver.function() << std::endl;
 
    //Selectively print the relaxed shape
    if(currPrintFlag){
      sstm << fname <<"-relaxed-" <<fileNum++;
      rName = sstm.str();
      
      model.print(rName);
      
      sstm.str("");
      sstm.clear(); // Clear state flags 
    }    

    // Calculating centre of the sphere as average of position vectors
    // of all nodes of the spherical capsid
    tvmet::Vector<double,3> Xavg(0.0);
    for ( int i = 0; i<defNodes.size(); i++){
      Xavg += defNodes[i]->point();
    }
    Xavg /= nodes.size();

    //We will calculate radius using the quadrature points
    std::vector<double> qpRadius;
    LSB::FeElementContainer fem = bd->shells();
    for(LSB::ConstFeElementIterator i=fem.begin(); i!=fem.end(); ++i){
      const LS::NodeContainer eleNodes = (*i)->nodes();      
      LS::QuadPointContainer qp = (*i)->quadraturePoints();
      for(LS::ConstQuadPointIterator qpi = qp.begin();qpi!=qp.end();++qpi){
	LoopShellShape s = (*qpi).shape;
	const LoopShellShape::FunctionArray fn = s.functions();
	tvmet::Vector<double,3> Xq(0.0);
	for(int i=0;i<fn.size();i++){
	  Xq += tvmet::mul(eleNodes[i]->point(),fn(i));
	}
	double qpR = tvmet::norm2(Xq-Xavg);
	qpRadius.push_back(qpR);
      }
    }
    
    Ravg = 0.0;
    for ( int i = 0; i<qpRadius.size(); i++){
       Ravg += qpRadius[i];
    }
    Ravg /= qpRadius.size();
  
    double dRavg2 = 0.0;
    for ( int i = 0; i<qpRadius.size(); i++){
      double dR =  qpRadius[i]-Ravg; 
      dRavg2 += dR*dR;
    }
    dRavg2 /= qpRadius.size();

    asphericity = dRavg2/(Ravg*Ravg);
    double gammaCalc = Y*Ravg*Ravg/KC;

    //Calculate Average Principal Strain
    std::vector<double> maxStrain =  bd->calcMaxPrincipalStrains();
    double avgStrain = 0.0;
    for(int e=0; e < maxStrain.size(); e++){
      avgStrain += maxStrain[e];
    }
    avgStrain /= maxStrain.size();

    myfile <<Ravg<<"\t"<<Y<<"\t"<<asphericity<<"\t"<<gammaCalc
	   <<"\t"<<avgStrain<<endl;

  }
  myfile.close();
  t2=clock();
  float diff ((float)t2-(float)t1);
  std::cout<<"Total execution time: "<<diff/CLOCKS_PER_SEC
	   <<" seconds"<<std::endl;
}
