#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include <tvmet/Vector.h>
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "C0MembraneBody.h"
#include "ShapeTri3.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "CgDescent.h"
#include "CgDescent-globals.h"
#include "ViscousRegularizer.h"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSetReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkWarpVector.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void computeAspherity(char * modelName, double & aspherity, double & radius);

int main(int argc, char* argv[])
{
  if( argc != 4 ) {
      cout << "Usage: relax modelName gamma_min gamma_max." << endl;
      return(0);
  }

  // Alert user to OpenMP status
#if defined(_OPENMP)
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;
#endif

  // Initialize MPI if in use
  bool verbose=true;
#ifdef WITH_MPI
  MPI_Init( &argc, &argv );
  int procId=0;
  MPI_Comm_rank( MPI_COMM_WORLD, &procId );
  if( procId !=0 ) verbose=false;
#endif

  ////////////////////////////////////////////////////////////////////
  // Input section
  ////////////////////////////////////////////////////////////////////

  // read in vtk file
  string modelName = argv[1];

  string inputFileName = modelName + ".vtk";

  vtkDataSetReader * reader = vtkDataSetReader::New();
  reader->SetFileName( inputFileName.c_str() );
  // send through normals filter to ensure that triangle orientations
  // are consistent
  vtkPolyDataNormals * normals = vtkPolyDataNormals::New();
  normals->SetInput( reader->GetOutput() );
  normals->ComputeCellNormalsOn();
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  vtkPolyData * mesh = normals->GetOutput();
  mesh->Update();
  std::cout << "mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
	    << std::endl;

  // create vector of nodes, get positions from vtk input
  double R = 1.0;
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
  if(verbose) cout << "Number of nodes: " <<nodes.size() << endl
		   << "Ravg = " << Ravg << endl;

  // read in connectivities
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

  // rescale reference such that average radius is R
  for(int i=0; i<defNodes.size(); i++) {
    DeformationNode<3>::Point x;
    x = defNodes[i]->position();
    x *= R/Ravg;
    defNodes[i]->setPosition(x);
    defNodes[i]->setPoint(x);
    mesh->GetPoints()->InsertPoint(i,&(x[0]));
  }



  ////////////////////////////////////////////////////////////////////
  // Solution section
  ////////////////////////////////////////////////////////////////////

  // open a file to output the FvK vs. aspherity data
  string aspName = modelName + ".aspherity.dat";
  ofstream asp(aspName.c_str());

  // loop over a range of FvK numbers, evenly spaced on the log axis
  double gammaMax = 1.0e10;
  double gammaMin = 10.0;
  gammaMin = atof( argv[2] );
  gammaMax = atof( argv[3] );
  std::cout << "Minimum FvK = " << gammaMin << std::endl;
  std::cout << "Maximum FvK = " << gammaMax << std::endl;
  for(double gamma=gammaMin; gamma <= 1.1*gammaMax; gamma *= sqrt(10.0)) {

    // initialize material response
    double Y = sqrt(gamma);
    double KC = 1.0/Y;
    double nu = 0.3;
    
    // Introduce a numerical scaling factor to make energy and forces
    // large enough to avoid issues with machine precision.  Later
    // divide forces by factor when printing to file.
    double scalingFactor=1.0e6;
    Y *= scalingFactor;
    KC *= scalingFactor;

    double KG = -KC;
    double C0 = 0.0;

    typedef FVK MaterialType;

    double Yb=Y;
    bool mixed=true;

    if(mixed) Yb=0.0;
    MaterialType bending( KC, KG, C0, Yb, nu );
   
    // create Body for bending
    unsigned int quadOrder = 2;
    typedef LoopShellBody<MaterialType> LSB;
    LSB bd(bending, connectivities, nodes, quadOrder);
    bd.setOutput(paraview);
    
    double viscosity = 0.0e-6;
    ViscousRegularizer vr(bd.nodes(), viscosity);
    bd.pushBack( &vr ); 
  
    double vrTol = 1.0e-10;

    // create Body for stretching
    std::vector< std::vector<int> > s_connectivities;
    for(int i=0; i<connectivities.size(); i++) {
      std::vector<int> c(3);
      for(int j=0; j<3; j++) c[j]=connectivities[i](j);
      s_connectivities.push_back(c);
    }

//     bd.resetReference(true);

    // create Model
    Model::BodyContainer bdc;
    bdc.push_back(&bd);

    MaterialType stretching( 0.0, 0.0, C0, Y, nu );
    quadOrder = 1;
    typedef C0MembraneBody<TriangleQuadrature,MaterialType,ShapeTri3> MB;
    MB bdm(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm.setOutput(paraview);

    // reset the ref configuration of all triangles to be equilateral
    bdm.resetEquilateral();

    if(mixed) {
      bdc.push_back(&bdm);
    }

    Model model(bdc,nodes);

//     // How long does a compute call take?
//     int numberOfComputes=100;
//     clock_t t1=clock();
//     for(int i=0; i<numberOfComputes; i++) bd.compute(true,true,false);
//     clock_t t2=clock();
//     std::cout << numberOfComputes << " compute calls took "
// 	      << (double)(t2-t1)/CLOCKS_PER_SEC
// 	      << " seconds."
// 	      << std::endl;
//     return 0;

//     // perturb current positions
//     srand(time(0));
//     for(int a=0; a<nodes.size(); a++) {
//       for(int i=0; i<nodes[a]->dof(); i++) {
// 	nodes[a]->addPoint(i,2e-2*((double)(rand())/RAND_MAX-0.5));
//       }
//     }
//     model.checkConsistency(true,false);
//     return 0;

#if 1
    // initialize Quasi-Newton BFGS solver
    int m=5;
    double factr=1.0e+7;
    double pgtol=1.0e-5;
    int iprint = 0;
    double pentol=1.0e-4;
    ifstream lbfgsbinp("lbfgsb.inp");
    lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
    if(verbose) 
      std::cout << "Input iprint: " << iprint << std::endl
		<< "Input factr: " << factr << std::endl
		<< "Input pgtol: " << pgtol << std::endl
		<< "Input m: " << m << std::endl
		<< "Input pentol: " << pentol << std::endl;
    Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );

    // set up bounds for solver
    bool boundaryConditionFlag = true;
    blitz::Array<int,1> nbd(3*nodes.size());
    blitz::Array<double,1> lo(3*nodes.size());
    blitz::Array<double,1> hi(3*nodes.size());
    nbd = 0;
    lo = 0.0;
    hi = 0.0;
    if(boundaryConditionFlag) {
      // constrain apices
      int a=0;
      for(int i=0; i<2; i++) {
	nbd(3*a+i) = 2;
	hi(3*a+i) = 0.0;
	lo(3*a+i) = 0.0;
      }
      a = 11;
      for(int i=0; i<2; i++) {
	nbd(3*a+i) = 2;
	hi(3*a+i) = 0.0;
	lo(3*a+i) = 0.0;
      }
    }


#else

    // set up CG solver
    double tol=1.0e-5;
    CgDescent solver(model.dof(), tol);

#endif

    // print starting energy, minimize, and output relaxed shape
    std::cout << "Relaxing shape for gamma = " << gamma << std::endl
	      << "Energy = " << solver.function() << std::endl;    


    double vrEnergy = vr.energy();
    double bdEnergy = bd.energy();
        int viter=0;
    do {
      if(verbose) std::cout << "viter = " << viter << std::endl;
      viter++;
      vr.step();
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bd.energy();
      if(verbose)
	std::cout << "v: " << vrEnergy << std::endl
		  << "e: " << bdEnergy << std::endl;
    } while(vrEnergy > vrTol * bdEnergy);

//     solver.solve( &model );

    std::cout << "Done relaxing shape." << std::endl;

    char fname[100];

    // update vtk mesh and print it

    vtkSmartPointer<vtkDoubleArray> displacements =
    vtkSmartPointer<vtkDoubleArray>::New();
    displacements->SetNumberOfComponents(3);
    displacements->SetNumberOfTuples(mesh->GetNumberOfPoints());
    displacements->SetName("Displacements");
 
    for(int a=0; a<mesh->GetNumberOfPoints(); a++) {      
      DeformationNode<3>::Point x(0.0);
      x = defNodes[a]->point() - defNodes[a]->position();
      displacements->SetTuple(a, &(x(0)) );
    }
    
    mesh->GetPointData()->AddArray(displacements);
 
    sprintf(fname,"%s-%g.vtk",modelName.c_str(), gamma);
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput( mesh );
    writer->SetFileName( fname ); 
    writer->SetFileTypeToBinary();
    writer->Write();
 
    sprintf(fname,"%s-%g-model",modelName.c_str(), gamma);
    model.print(fname);
   
    // calculate new average radius and aspherity from control points
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
    double aspherity = dRavg2/(Ravg*Ravg);

    // Figure out whether elements are colored by valence (penton/hexon)
    vtkCellData *celldata = mesh->GetCellData();
    if (!celldata) {
      cout << "No cell data in input mesh." << endl;
    }

    std::cout << " contains cell data with " 
	      << celldata->GetNumberOfArrays() 
	      << " arrays." << std::endl;
    std::string arrayname;
    for (int i = 0; i < celldata->GetNumberOfArrays(); i++) {
      arrayname = 
	(celldata->GetArrayName(i) ? celldata->GetArrayName(i) : "NULL") ;
      if( arrayname == "valence" || arrayname == "eletype" ) {
	std::cout << "\tArray " << i 
		  << " is named "
		  << arrayname
		  << std::endl;
	break;
      } else {
	cout << "Couldn't find valence or eletype cellData array in mesh." 
	     << endl;
      }
    }

    //Threshold out pentons
    vtkSmartPointer<vtkThreshold> threshold = 
      vtkSmartPointer<vtkThreshold>::New();
    threshold->SetInput(mesh);
    threshold->ThresholdByUpper(5.5);
    threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayname.c_str());
    threshold->Update();
 
    vtkUnstructuredGrid* thresholdedPolydata = threshold->GetOutput();
    std::cout << "There are " << thresholdedPolydata->GetNumberOfCells() 
	      << " cells after thresholding." << std::endl;

    thresholdedPolydata->GetPointData()->SetActiveVectors("Displacements");

    vtkSmartPointer<vtkWarpVector> warp = vtkSmartPointer<vtkWarpVector>::New();
    warp->SetInput(thresholdedPolydata);
    warp->SetScaleFactor(1.0);

    warp->Update();

    vtkPointSet * warpedGrid = warp->GetOutput();


    // sprintf(fname,"%s-%g-thresh.vtk",modelName.c_str(), gamma);
    // vtkUnstructuredGridWriter *gridwriter = vtkUnstructuredGridWriter::New();
    // gridwriter->SetInput( warp->GetOutput() );
    // gridwriter->SetFileName( fname ); 
    // gridwriter->SetFileTypeToBinary();
    // gridwriter->Write();

    Vector3D XavgThresh(0.0);
    for ( int i = 0; i<warpedGrid->GetNumberOfPoints(); i++){
      Vector3D Xi( warpedGrid->GetPoint(i), warpedGrid->GetPoint(i)+3 );
      XavgThresh += Xi;
    }
    XavgThresh /= nodes.size();
    double RavgThresh = 0.0;
    for ( int i = 0; i<warpedGrid->GetNumberOfPoints(); i++){
      Vector3D Xi( warpedGrid->GetPoint(i), warpedGrid->GetPoint(i)+3 );
      RavgThresh += tvmet::norm2( Xi - XavgThresh );
    }
    RavgThresh /= warpedGrid->GetNumberOfPoints();
    
    double dRavg2Thresh = 0.0;
    for ( int i = 0; i<warpedGrid->GetNumberOfPoints(); i++){
      Vector3D Xi( warpedGrid->GetPoint(i), warpedGrid->GetPoint(i)+3 );
      double dR = (tvmet::norm2( Xi - XavgThresh) - RavgThresh); 
      dRavg2Thresh += dR*dR;
    }
    dRavg2Thresh /= warpedGrid->GetNumberOfPoints();

    double gammaThresh = Y*RavgThresh*RavgThresh/KC;
    double aspherityThresh = dRavg2Thresh/(RavgThresh*RavgThresh);

    // Estimate aspherity of limit surface by subdividing a few times
    // with VTK
//     double aspherity_lim = 0.0;
//     double radius_lim = 0.0;
//     computeAspherity( fname, aspherity_lim, radius_lim );
//     double gamma_lim = Y*radius_lim*radius_lim/KC;
    double Rarea = sqrt( bd.area()/(4.0*M_PI) );
    if(verbose) 
      cout << "            Ravg = " << Ravg << endl
	   << "      RavgThresh = " << RavgThresh << endl
	   << "           Rarea = " << Rarea << endl
	   << "            Xavg = " << Xavg << endl
	   << "      XavgThresh = " << XavgThresh << endl
	   << "      FvK number = " << gammaCalc << endl
	   << "      FvK Thresh = " << gammaThresh << endl
	   << "       Aspherity = " << aspherity << endl
	   << "Aspherity Thresh = " << aspherityThresh << endl;
// 	   << "Ravg(limit) = " << radius_lim << endl
// 	   << "FvK(limit) = " << gamma_lim << endl
// 	   << "Aspherity(limit) = " << aspherity_lim << endl;
    asp << std::setw(18) << gammaCalc
	<< std::setw(18) << aspherity
	<< std::setw(18) << gammaThresh
	<< std::setw(18) << aspherityThresh
	<< std::endl; 
  }
  
  std::cout << "All done.  Bye now." << std::endl;
  return 0;
}


//
// Estimate aspherity of limit surface by subdividing a few times with VTK
//
void computeAspherity(char * modelName, double & aspherity, double & radius ) {

  string fileName = string(modelName) + "-bd0.vtk";
  
  aspherity=0;

  vtkDataSetReader *reader = vtkDataSetReader::New();
  reader->SetFileName( fileName.c_str() );
  reader->SetVectorsName("displacements");

  vtkWarpVector * warp = vtkWarpVector::New();
  warp->SetInput( reader->GetOutput() );
  warp->SetScaleFactor(1.0);

  int n = 1;
  vtkLoopSubdivisionFilter *subdivisionFilter 
    = vtkLoopSubdivisionFilter::New();
  subdivisionFilter->SetNumberOfSubdivisions( n );
  subdivisionFilter->SetInput( warp->GetOutput() );
  subdivisionFilter->Update();

//   vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
//   writer->SetInput( subdivisionFilter->GetOutput() );
//   writer->SetFileName( "junk.vtk" ); 
//   writer->SetFileTypeToASCII();
//   writer->Write();
  
//   ifstream ifs( "junk.vtk" );
//   if (!ifs) {
//     cout << "Cannot open input file: junk.vtk" << endl;
//     exit(0);
//   }

  // find points header
  int npts= warp->GetOutput()->GetNumberOfPoints();
  std::cout << "Number of points before subdivision: " << npts << std::endl;
  vtkPolyData * mesh = subdivisionFilter->GetOutput();
  npts=mesh->GetNumberOfPoints();
  std::cout << "Number of points after subdivision:  " << npts << std::endl;

  std::vector<Vector3D> points;
  points.reserve(npts);

  // read in points
  for(int a=0; a<npts; a++) {
    Vector3D x;
    mesh->GetPoint(a, &(x[0]));
    points.push_back( x );
  }

  subdivisionFilter->Delete();
  warp->Delete();
  reader->Delete();

//   while( token != "displacements" ) ifs >> token;
//   ifs >> token;// skip number type
//   for(int i=0; i<npts; i++) {
//     Vector3D u;
//     ifs >> u(0) >> u(1) >> u(2);
//     points[i]+=u;
//   }

  Vector3D Xavg(0.0);
  for ( int i = 0; i<npts; i++){
    Xavg += points[i];
  }
  Xavg /= npts;
    for ( int i = 0; i<npts; i++){
    radius += tvmet::norm2( points[i] - Xavg );
  }
  radius /= npts;
    
  double dRavg2 = 0.0;
  for ( int i = 0; i<npts; i++){
    double dR = (tvmet::norm2( points[i] - Xavg) - radius); 
    dRavg2 += dR*dR;
  }
  dRavg2 /= npts;

  aspherity = dRavg2/(radius*radius);

  return;
}
