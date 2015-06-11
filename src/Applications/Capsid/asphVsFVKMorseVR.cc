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
#include "Quadrature.h"
#include "TriangleQuadrature.h"

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
#include "PotentialBody.h"
#include "Utils/PrintingProtein.h"
#include "ViscousRegularizer.h"

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
  string inputFileName = argv[1];

  bool remesh = true;

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
  double Ravg = 0.0;

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

  // Rescale size of the capsid
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = defNodes[i]->point();
    x *= Rcapsid/Ravg;
    defNodes[i]->setPoint(x);
    defNodes[i]->setPosition(x);
  }
  Ravg = 1.0;//At least, we expect it to be 1.0 now
  
  // Calculate side lengths of the equilateral triangles
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
  std::cout<<"Average equilateral triangle edge length:"
	   <<EquilateralEdgeLength<<endl;

  //******************* READ FVK DATA FROM FILE *************************//

  
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
  
  //Uncomment the block below to calculate FVK steps by code
  /*
  // FVK = 10^(exponent)
  double expMin = 0.0; // minimum exponent -> FVK = 10^0
  double expMax = 4.0; // maximum exponent -> FVK = 10^4
  double FVK_steps = 1000;
  double expIncr = (expMax - expMin)/FVK_steps;
  double currExponent;
  */

  std::stringstream sstm;
  string fname = inputFileName.substr(0,inputFileName.find("."));
  string iName;
  string rName;

  ofstream myfile;
  myfile.open ("asphVsFVKMorse.dat");
  myfile << "Ravg,Y,asphericity,FVKin,FVKout" << endl;

  
  //Parameters for the l-BFGS solver
  int m=5;
  int maxIter=1e6;
  double factr=1.0e+1;
  double pgtol=1e-7;
  int iprint = 1;
  Lbfgsb solver(3*nodes.size(), m, factr, pgtol, iprint, maxIter );
  
  typedef FVK MaterialType;
  typedef LoopShellBody<MaterialType> LSB;
  typedef LoopShell<MaterialType> LS;

  //***************************  SOLUTION LOOP ***************************//

  //Loop over all values of gamma and relax the shapes to get
  //asphericity

  for(std::vector<vector<double> >::iterator q=gammaVec.begin();
      q!=gammaVec.end();++q){

    /*
      currExponent = expMin + q*expIncr;
      double gamma = pow(10.0,currExponent);
    */

    double gamma = (*q)[0];
    double currPrintFlag = (*q)[1];    
    
    //****************  Protein body parameters ****************//

    //For Morse material 
    double Rshift = EquilateralEdgeLength;
    double ratio = 5.0;
    double epsilon = 0.0023/ratio;
    double sigma = (6.0/Rshift)*sqrt(ratio);
    double springConstant = 2*sigma*sigma*epsilon;
    double ARtol = 1.5;
    double PotentialSearchRF=1.5*Ravg;    

    //*************** Bending body parameters ***************//
    double Y = 2.0/sqrt(3)*springConstant; // Young's modulus
    double nu = 1.0/3.0;
    double KC = Y*Ravg*Ravg/gamma; // Bending modulus calculated from gamma and Y
    double KG = -2*(1-nu)*KC; // Gaussian modulus
    double C0 = 0.0;
    int quadOrder = 2;
    
    //**************** The Bodies *****************//

    MaterialType bending(KC,KG,C0,0.0,0.0);
    LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
    bd->setOutput(paraview);

    Morse Mat(epsilon,sigma,Rshift);
    PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);

    // add some viscosity for regularization
    double viscosity_inp = 20;
    double minViscosity = 1.0e-6*viscosity_inp;
    double maxViscosity = 1.0e+6*viscosity_inp;
    ViscousRegularizer vr(bd->nodes(), viscosity_inp);
    bd->pushBack( &vr );

    // set viscosity parameters
    double targetVelocity = 1e-4*Ravg;
    double vrTol = 1.0e-10;

    double vrEnergy = vr.energy();
    double bdEnergy = bd->energy();
    double PrBodyEnergy = PrBody->energy();

    //********************** Create Model ********************//
    Model::BodyContainer bdc;
    bdc.push_back(PrBody);
    bdc.push_back(bd);    
    Model model(bdc,nodes);
     
    std::cout<< "Morse potential parameters:" << endl
	     << "sigma = " << sigma <<" epsilon = " << epsilon
	     << " Rshift = "<< Rshift <<endl;
     
    PrBody->compute(true, false, false);
    std::cout << "Initial protein body energy = " << PrBody->energy() << endl;
     
    std::cout << "Relaxing shape for gamma = " << gamma<< std::endl
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

    int viterMax = 20;
    for(int viter=0; viter<viterMax; viter++) {

      if(viter==viterMax-1) vr.setViscosity(minViscosity);

      blitz::Array<double,1> vSave(model.dof());
      model.getField(solver);
      for(int i=0; i<model.dof(); i++ ) vSave(i) = solver.field(i);
      
      std::cout << std::endl 
		<< "VISCOUS ITERATION: " << viter 
		<< "\t viscosity = " << vr.viscosity()
		<< std::endl
		<< std::endl;
  
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bd->energy();
      PrBodyEnergy = PrBody->energy();

      std::cout << "ENERGY:" << std::endl
		<< "viscous energy = " << vrEnergy << std::endl
		<< " LSB  body energy = " << bdEnergy << std::endl
		<< " Protein  body energy = " << PrBodyEnergy << std::endl
		<< "  total energy = " << solver.function() << std::endl
		<< std::endl;
	
      std::cout << "VISCOSITY: " << std::endl
		<< "          velocity = " << vr.velocity() << std::endl
		<< "   target velocity = " << targetVelocity << std::endl
		<< " updated viscosity = " << vr.viscosity() << std::endl
		<< std::endl;

      // step forward in "time", relaxing viscous energy & forces 
      vr.step();

      if(vrEnergy < std::abs(vrTol*bdEnergy) && solver.projectedGradientNorm()<=pgtol) {
	// viscous energy is small enough; exit
	break;
      }
    }     

    //Uncomment the following region if you want to print initial
    //shapes
    /*
      sstm << fname <<".initialFVK-" <<gamma;
      iName = sstm.str();
       
      model.print(iName);
       
      sstm.str("");
      sstm.clear(); // Clear state flags.
    */
     
    std::cout << "Shape relaxed." << std::endl
	      << "Energy = " << solver.function() << std::endl;
     
    //Selectively print the relaxed shapes
    if(currPrintFlag){
      sstm << fname <<".relaxedFVK-" << gamma;
      rName = sstm.str();

      model.print(rName);

      sstm.str("");
      sstm.clear(); // Clear state flags    
    }

    //Re-calculate average equilateral triangle edge length
    double TriangleEdgeLength = 0.0;
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
      TriangleEdgeLength += 
	(tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))/3.0;
    }
    TriangleEdgeLength /= connectivities.size();
    std::cout<<"Average equilateral triangle edge length after relaxing:"
	     <<TriangleEdgeLength<<endl;
    
    tvmet::Vector<double,3> Xavg(0.0);
    for ( int i = 0; i<defNodes.size(); i++){
      Xavg += defNodes[i]->point();
    }
    Xavg /= defNodes.size();
    
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

    double asphericity = dRavg2/(Ravg*Ravg);
    double gammaCalc = Y*Ravg*Ravg/KC;
  
    myfile<<Ravg<<","<<Y<<","<<asphericity
	  <<","<<gamma<<","<<gammaCalc<<endl;
    
    //Release the dynamically allocated memory
    delete bd;
    delete PrBody;
  }

  myfile.close();
  t2=clock();
  float diff ((float)t2-(float)t1);
  std::cout<<"Total execution time: "<<diff/CLOCKS_PER_SEC
	   <<" seconds"<<std::endl;
}
    
