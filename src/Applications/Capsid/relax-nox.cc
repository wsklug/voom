#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include <tvmet/Vector.h>

#include "NOX.H"
#include "NOX_LAPACK_Group.H"

#include <vtk/vtkCell.h>
#include <vtk/vtkDataSetReader.h>
#include <vtk/vtkPolyDataWriter.h>
#include <vtk/vtkPolyDataNormals.h>
#include <vtk/vtkWarpVector.h>
#include <vtk/vtkLoopSubdivisionFilter.h>

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



#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

void computeAspherity(char * modelName, double & aspherity, double & radius);

class Voom2NOX : public NOX::LAPACK::Interface {

public:

  typedef std::vector< Body* > 			BodyContainer;
  typedef std::vector< Body* >::iterator 	BodyIterator;
  typedef std::vector< Body* >::const_iterator 	ConstBodyIterator;
  
  typedef std::vector< NodeBase* > 			NodeContainer;
  typedef std::vector< NodeBase* >::iterator 		NodeIterator;
  typedef std::vector< NodeBase* >::const_iterator 	ConstNodeIterator;
  
  //! Constructor
  Voom2NOX(BodyContainer & bodies, NodeContainer & nodes, int nDOF) 
    : _bodies(bodies), _nodes(nodes) , _initialGuess(nDOF), _solution(nDOF)
  { }
    
  //! Destructor
  ~Voom2NOX() {};

  const NOX::LAPACK::Vector& getInitialGuess()
  {
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      const NodeBase::DofIndexMap & idx = (*n)->index();
      for(int ni=0; ni<(*n)->dof(); ni++)
	_initialGuess(idx[ni]) = (*n)->getPoint(ni);
    }
    return _initialGuess;
  };

  //! Return true solution vector
  const NOX::LAPACK::Vector& getSolution()
  {
    return _solution;
  };

  bool computeF(NOX::LAPACK::Vector& f, const NOX::LAPACK::Vector &x)
  {
    // zero out all forces and stiffness in nodes before computing bodies
    
    for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      for(int i=0; i<(*n)->dof(); i++) (*n)->setForce(i,0.0);
    }

    // send trial solution to bodies
    for(NodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      const NodeBase::DofIndexMap & idx = (*n)->index();
      for(int ni=0; ni<(*n)->dof(); ni++)
	(*n)->setPoint(ni, x(idx[ni]));
    }

    bool f0=false;
    bool f1=true;
    bool f2=false;

    // compute bodies
    for(BodyIterator b=_bodies.begin(); b!=_bodies.end(); b++)
      (*b)->compute( f0, f1, f2);

    // zero out force array
    f.init(0.0);

    // assemble forces
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) {
      const NodeBase::DofIndexMap & idx = (*n)->index();
	for(int i=0; i<idx.size(); i++)
	  f( idx[i] ) = (*n)->getForce(i);
    }

    return true;
  };
  
  bool computeJacobian(NOX::LAPACK::Matrix<double>& J, 
		       const NOX::LAPACK::Vector & x)
  {
    throw "No Jacobian.";
    return false;
  };

private:

  //! Container of bodies that will compute forces
  BodyContainer _bodies;

  //! Container of nodes where forces will be assembled
  NodeContainer _nodes;

  //! Initial guess
  NOX::LAPACK::Vector _initialGuess;
  //! Correct solution
  NOX::LAPACK::Vector _solution;

};


int main(int argc, char* argv[])
{
  if( argc != 2 ) {
      cout << "Usage: relax modelName." << endl;
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
  normals->ConsistencyOn();
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
//   normals->Update();
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
//     x = defNodes[i]->point();
//     x *= R/tvmet::norm2(x);
    defNodes[i]->setPoint(x);
  }



  ////////////////////////////////////////////////////////////////////
  // Solution section
  ////////////////////////////////////////////////////////////////////

  // open a file to output the FvK vs. aspherity data
  string aspName = modelName + ".aspherity.dat";
  ofstream asp(aspName.c_str());

  // loop over a range of FvK numbers, evenly spaced on the log axis
  double gammaMax = 1.0e3;
  double gammaMin = 1.0e3;
  for(double gamma=gammaMax; gamma >= gammaMin; gamma /= sqrt(10.0)) {

    // initialize material response
    double Y = sqrt(gamma);
    double KC = 1.0/Y;
    double nu = 0.3;
    
    double KG = -KC;
    double C0 = 0.0;

    typedef FVK MaterialType;

    double Yb=Y;
    bool mixed=false;

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
    Voom2NOX::BodyContainer bdc;
    bdc.push_back(&bd);

    MaterialType stretching( 0.0, 0.0, C0, Y, nu );
    quadOrder = 1;
    typedef C0MembraneBody<TriangleQuadrature,MaterialType,ShapeTri3> MB;
    MB bdm(stretching, s_connectivities, nodes, quadOrder, 0.0);
    bdm.setOutput(paraview);
    if(mixed) {
      bdc.push_back(&bdm);
    }

    // Set up NOX solver interface
    Voom2NOX noxInterface(bdc,nodes,dof);

    // Create a group which uses that problem interface. The group will
    // be initialized to contain the default initial guess for the
    // specified problem.
    Teuchos::RCP<NOX::LAPACK::Group> grp = 
      Teuchos::rcp(new NOX::LAPACK::Group(noxInterface));

    // Set up the status tests
    Teuchos::RCP<NOX::StatusTest::NormF> statusTestA = 
      Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-4));
    Teuchos::RCP<NOX::StatusTest::MaxIters> statusTestB = 
      Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
    Teuchos::RCP<NOX::StatusTest::Combo> statusTestsCombo = 
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
					      statusTestA, statusTestB));
    
    // Create the list of solver parameters
    Teuchos::RCP<Teuchos::ParameterList> solverParametersPtr =
      Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& solverParameters = *solverParametersPtr;
    
    // Set the level of output (this is the default)
    solverParameters.sublist("Printing").set("Output Information", 
					     NOX::Utils::Warning + 
					     NOX::Utils::OuterIteration + 
					     NOX::Utils::InnerIteration + 
					     NOX::Utils::Parameters);
    
    // Set the solver (this is the default)
    solverParameters.set("Nonlinear Solver", "Line Search Based");

    // Create the direction parameters sublist
    Teuchos::ParameterList& directionParameters = solverParameters.sublist("Direction");
    
    directionParameters.set("Method","Steepest Descent");
    
    // Create the line search parameters sublist
    Teuchos::ParameterList& lineSearchParameters = solverParameters.sublist("Line Search");
    
    // Set the line search method
    lineSearchParameters.set("Method","More'-Thuente");
    
    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = 
      NOX::Solver::buildSolver(grp, statusTestsCombo, solverParametersPtr);
    
    // Solve the nonlinesar system
    NOX::StatusTest::StatusType status = solver->solve();
    
    // Warn user if solve failed
    if (status == NOX::StatusTest::Converged)
      cout << "Solver converged!" << endl;
    else
      cout << "Error: Solve failed to converge!" << endl;
    
    // Print the parameter list
    cout << "\n" << "-- Parameter List From Solver --" << "\n";
    solver->getList().print(cout);
    
    
    // Get the answer
    NOX::LAPACK::Group solnGrp = 
      dynamic_cast<const NOX::LAPACK::Group&>(solver->getSolutionGroup());
    
    // Print the answer
//     cout << "\n" << "-- Final Solution From Solver --" << "\n";
//     solnGrp.print();
    
    // Print the expected answer
//     solnGrp.setX(noxInterface.getSolution());
//     solnGrp.computeF();
    
    char fname[100];
    sprintf(fname,"%s-%g",modelName.c_str(), gamma);
    bdc[0]->printParaview(fname);

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

    // Estimate aspherity of limit surface by subdividing a few times
    // with VTK
    double aspherity_lim = 0.0;
    double radius_lim = 0.0;
    computeAspherity( fname, aspherity_lim, radius_lim );
    double gamma_lim = Y*radius_lim*radius_lim/KC;
    double Rarea = sqrt( bd.area()/(4.0*M_PI) );
    if(verbose) 
      cout << "Ravg = " << Ravg << endl
	   << "Rarea = " << Rarea << endl
	   << "Ravg(limit) = " << radius_lim << endl
	   << "Xavg = " << Xavg << endl
	   << "FvK number = " << gammaCalc << endl
	   << "FvK(limit) = " << gamma_lim << endl
	   << "Aspherity = " << aspherity << endl
	   << "Aspherity(limit) = " << aspherity_lim << endl;
    asp << std::setw(18) << gamma 
	<< std::setw(18) << aspherity
	<< std::setw(18) << gamma_lim
	<< std::setw(18) << aspherity_lim << std::endl; 
  }
  
  std::cout << "All done.  Bye now." << std::endl;
  return 0;
}


//
// Estimate aspherity of limit surface by subdividing a few times with VTK
//
void computeAspherity(char * modelName, double & aspherity, double & radius ) {

  string fileName = string(modelName) + ".vtk";
  
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
