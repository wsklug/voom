#include <string>
#include <iostream>
#include <vector>
#include <fstream>
//#include <getopt.h>
#include "Node.h"
#include "StVenant.h"
#include <tvmet/Vector.h>
#include "FiniteContractionBody.h"
#include "ContractionWrapper.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include "Model.h"
#include "Solver.h"
#include "Lbfgsb.h"
#include "SimulatedAnnealing.h"
//#include "Rigid.h" 

using namespace tvmet;
using namespace std;
using namespace voom;

#if defined(_OPENMP)
#include <omp.h>
#endif

int main(int argc, char* argv[])
{

#ifdef _OPENMP
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


    bool boundaryConditionFlag = false;

  //
  // create vector of nodes
  int dof=0;
  
  int npts=4;
  std::vector< NodeBase* > nodes;
  //std::vector< NodeBase* > defNodes;
  std::vector< DeformationNode<3>* > defNodes;
  
  nodes.reserve(npts);
  defNodes.reserve(npts);
  
  std::vector< ScalarFieldNode<3>* > voltNodes;
  voltNodes.reserve(npts);

  double volts = 0.0;
    
  std::vector<Vector3D> X(4);
  X[0] = 0.0, 0.0, 0.0;
  X[1] = 2.0, 0.0, 0.0;
  X[2] = 0.0, 2.0, 0.0;
  X[3] = 0.0, 0.0, 2.0;

  srand(time(0));

  for(int id=0; id<4; id++) {

    Vector3D x; // current (deformed) position
    x = X[id];

    // perturb current positions by random displacement
    for(int i=0; i<3; i++) {
//       x(i) += 0.1*rand()/((double)RAND_MAX);
      x(i) *= 1.5;
    }
    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    defNodes.push_back( new DeformationNode<3>(id,idx,X[id],x) );
    // have to put deformation nodes into nodes[] array so that solver
    // can minimize energy with respect to their current positions.
    nodes.push_back( defNodes[id] );

    voltNodes.push_back( new ScalarFieldNode<3>(id,idx,X[id],volts) );
    
  }

// connectivities  
 
//    std::vector< std::vector<int,4> > connectivities;
//	std::vector<int, 4> c;
     std::vector< std::vector<int> > connectivities;
       std::vector<int> c(4);
	connectivities.reserve(1);
	c[0] = 1;
	c[1] = 2;
	c[2] = 3;
	c[3] = 4;
	connectivities.push_back(c);
 

  // create material
  double rho = 500.0;
  double E = 200.0;//1.0e3;
  double nu = 0.4;
  typedef StVenant MaterialType;
  MaterialType muscle( rho, E, nu );
    
  typedef ContractionWrapper<StVenant> ContractionMaterial;
  ContractionMaterial cardiacMaterial(muscle);
  std::cout << "Material has been created." << std::endl;
	
  // create Body
  const unsigned quadOrder = 2;
  
  typedef FiniteContractionBody<TetQuadrature, ContractionMaterial, ShapeTet4,
  TetQuadrature, ContractionMaterial, ShapeTet4> BodyType; 

  BodyType bd(cardiacMaterial, cardiacMaterial, connectivities, defNodes, voltNodes, quadOrder);




  //bd.setOutput(Body::paraview);

  bd.compute(true,true,false);
//  bd.printParaviewQuadTet("initial");
  std::cout << "body energy = " << bd.energy() << std::endl;

  //return 0;
  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  Model model(bdc, nodes);

//   model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
//   return 0;

  int mm=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = -1;
//  ifstream lbfgsbinp("lbfgsb.inp");
//  lbfgsbinp >> iprint >> factr >> pgtol >> mm ;
//  if(verbose) 
//    std::cout << "Input iprint: " << iprint << std::endl
//	      << "Input factr: " << factr << std::endl
//	      << "Input pgtol: " << pgtol << std::endl
//	      << "Input m: " << mm << std::endl;
  std::cout << "model.dof()"<<model.dof()<< std::endl;
  Lbfgsb solver(model.dof(), mm, factr, pgtol, iprint );//(true);

  
  // set up bounds for solver
  blitz::Array<int,1> nbd(3*nodes.size());
  blitz::Array<double,1> lo(3*nodes.size());
  blitz::Array<double,1> hi(3*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;

  solver.setBounds(nbd,lo,hi); 

  solver.solve( &model );

  std::cout << "All done.  Bye now." << std::endl;

  return 0;
}
