#include <string>
#include <iostream>
#include <time.h>
#include <vector>
#include <fstream>

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
  cout << "Number of nodes: " <<nodes.size() << endl;

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
  Ravg = 1.0;

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
  
  // Bending body material and other properties
  double Y;
  double KC;
  double KG;
  double C0 = 0.0; 
  int quadOrder = 2;
  typedef FVK MaterialType;
  typedef LoopShellBody<MaterialType> LSB;  
  KC = 1.0;
  KG = -2.0*(1.0-(1.0/3.0))*KC;
  
  // We will set Yb and nu as 0 so that LSB does not do stretching
  MaterialType bending(KC,KG,C0,0.0,0.0);
  LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
  bd->setOutput(paraview);

  // Protein body properties
  double PotentialSearchRF=1.5;  
  //For Lennard-Jones sigma = a/(2^(1/6)) where a = EquilateralEdge
  double sigma = EquilateralEdgeLength/1.122462048;
  double epsilon;
  double ARtol = 1.5; 

  // Initiliaze protein material
  LennardJones Mat(0.0,sigma); // We will update epsilon in the for-loop

  // Then initialize potential body
  PotentialBody * PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);    

  //Create the model
  Model::BodyContainer bdc;
  bdc.push_back(bd);
  bdc.push_back(PrBody);
  Model model(bdc,nodes);

  //Parameters for the l-BFGS solver
  int m=5;
  int maxIter=1000;
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = 1;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );//(true);

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

  //Uncomment the block comment to calculate FVK by code instead of
  //reading from dat file
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
  myfile.open ("asphVsFVK.dat");
  myfile << "Epsilon,Sigma,Ravg,Y,asphericity,FVK"<< endl;

  double asphericity;
  double gammaCalc;

  //***************************  FOR_LOOP ***************************//

  //Loop over all values of gamma and relax the shapes to get
  //asphericity
  /*for(int q=0;q <= FVK_steps; q++){
      currExponent = expMin + q*expIncr;
      double gamma = pow(10.0,currExponent);*/
  for(std::vector<vector<double> >::iterator q=gammaVec.begin();
      q!=gammaVec.end();++q){  
  
    double gamma = (*q)[0];
    double currPrintFlag = (*q)[1];

    //Update the epsilon for material
    Y = gamma/(Ravg*Ravg);
    epsilon = 0.0303089898881*sigma*sigma*Y;
    Mat.setEpsilon(epsilon);
  
    std::cout<< "Lennard-Jones potential parameters:" << endl
	     << "sigma =" << sigma <<" epsilon =" << epsilon << endl;
  
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
  
    //Uncomment the following region if you want to print initial
    //shapes
    /*
      sstm << fname <<".initial-" << q;
      iName = sstm.str();
    
      model.print(iName);
    
      sstm.str("");
      sstm.clear(); // Clear state flags.
    */
  
    // relax initial shape
    solver.solve(&model);
  
    std::cout << "Shape relaxed." << std::endl
	      << "Energy = " << solver.function() << std::endl;
 
    //Print relaxed shapes for selected FVK
    if(currPrintFlag){
      sstm << fname <<".relaxedFVK-" << gamma;
      rName = sstm.str();

      model.print(rName);

      sstm.str("");
      sstm.clear(); // Clear state flags    
    }

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

    gammaCalc = Y*Ravg*Ravg/KC;
    asphericity = dRavg2/(Ravg*Ravg);
  
    myfile<<epsilon<<","<<sigma<<","<<Ravg<<","<<Y<<","
	  <<asphericity<<","<<gammaCalc<<endl;
  }
  myfile.close();
  t2=clock();
  float diff ((float)t2-(float)t1);
  std::cout<<"Total execution time: "<<diff/CLOCKS_PER_SEC
	   <<" seconds"<<std::endl;
}
    
