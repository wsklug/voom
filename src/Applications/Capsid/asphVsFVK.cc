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

  double gamma_inp=1500.0;
  bool remesh = true;
  string prestressFlag = "yes";

  string inputFileName = "T7input.vtk";
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

  double Y = sqrt(gamma_inp);
  double KC = 1.0/Y;

  // Introduce a numerical scaling factor to make energy and forces
  // large enough to avoid issues with machine precision.  Later
  // divide forces by factor when printing to file.
  double scalingFactor = 1.0;//e6;
  Y *= scalingFactor;
  KC *= scalingFactor;

  //Amit: Set nu=0 and Yb=0 so that LoopShellBody does not handle
  //stretching anymore.
  //double nu = 1.0/3.0;
  //double Yb=Y;
  double Yb=0.0;  
  double nu = 0.0;
  double KG = -2.0*(1.0-nu)*KC;

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

  //Parameters for the l-BFGS solver
  int m=5;
  int maxIter=500;//1000;//model.dof();
  double factr=1.0e+1;
  double pgtol=1.0e-5;
  int iprint = 1;

  // Potential input parameters
  double PotentialSearchRF;  
  double sigma;
  double ARtol = 1.5;

  //For Lennard-Jones sigma = a/(2^(1/6)) where a = EquilateralEdgeLength  
  sigma = EquilateralEdgeLength/1.122462048;
  
  //We will look for only neighbouring nodes 
  PotentialSearchRF = 1.5*EquilateralEdgeLength;

  double FVK_min = 1.0;
  double FVK_max = 10000.0;
  double epsilon_min = 0.03030898988810338*(sigma/Ravg)*(sigma/Ravg)*KC*FVK_min;
  double epsilon_max = 0.03030898988810338*(sigma/Ravg)*(sigma/Ravg)*KC*FVK_max;
  double d_epsilon = (epsilon_max - epsilon_min)/2000;
  double epsilon;
 
  std::stringstream sstm;
  string fname = "T7input";
  string iName;
  string rName;

  ofstream myfile;
  myfile.open ("asphVsFVK.dat");
  myfile << "Epsilon,Sigma,Ravg,Y,asphericity,FVK"<< endl;

  double Ycalc;
  double asphericity;
  double gammaCalc;

  PotentialBody * PrBody;

  // Protein body implemented using Lennard-Jones body
  // Initiliaze potential material
  LennardJones Mat(epsilon_min, sigma);
  
  // Then initialize potential body
  PrBody = new PotentialBody(&Mat, defNodes, PotentialSearchRF);  
   
  bdc.push_back(PrBody);   
  Model model(bdc,nodes);  
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxIter );//(true);
  

  //***************************  FOR_LOOP ***************************//

  //Loop over all values of epsilon and relax the shapes to get
  //asphericity and FVK numbers
  for(int q=0;q <= 2000; q++){    
    epsilon = epsilon_min + q*d_epsilon;

    //Update the epsilon for material
    Mat.setEpsilon(epsilon);

    std::cout<< "Lennard-Jones potential parameters:" << endl
	     << "sigma =" << sigma <<" epsilon =" << epsilon << endl;
	      
    PrBody->compute(true, false, false);
    std::cout << "Initial protein body energy = " << PrBody->energy() << endl;
    
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
    
    //Uncomment the following region if you want to print relaxed
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
 
    //Uncomment the following region if you want to print relaxed
    //shapes
    /*
    sstm << fname <<".relaxed-" << q;
    rName = sstm.str();

    model.print(rName);

    sstm.str("");
    sstm.clear(); // Clear state flags    
    */

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

    Ycalc =  32.993511288*epsilon/(sigma*sigma);
    gammaCalc = Ycalc*Ravg*Ravg/KC;
    asphericity = dRavg2/(Ravg*Ravg);
  
    myfile<<epsilon<<","<<sigma<<","<<Ravg<<","<<Ycalc<<","
	  <<asphericity<<","<<gammaCalc<<endl;
  }
  myfile.close();
}
    
