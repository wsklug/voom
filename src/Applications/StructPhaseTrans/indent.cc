#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include "GLElastic.h"
#include "GLBody.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Rigid.h"
#include "ContactSPT.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

class ViscousRegularizer : public Element {

public:
  using Element::BaseNodeContainer; 
  using Element::BaseNodeIterator; 
  using Element::ConstBaseNodeIterator; 

  ViscousRegularizer( const BaseNodeContainer & nodes, double viscosity ) 
  {
    _baseNodes = nodes;
    _viscosity = viscosity;
    _energy = 0.0;
    int dof=0;
    for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      dof+=(*n)->dof();
    }
    _reference.resize(dof);
    _reference = 0.0;
    step();
  }

  void step() {
    int I=0;
    for(ConstBaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      for(int i=0; i<(*n)->dof(); i++,I++) {
	_reference(I) = (*n)->getPoint(i);
      }
    }
  }

  void compute(bool f0, bool f1, bool f2) {
    if(f0) _energy = 0.0;
    int I=0;
    for(BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
      for(int i=0; i<(*n)->dof(); i++,I++) {
	double dx = (*n)->getPoint(i) - _reference(I);
	    if(f0) _energy += 0.5 * _viscosity * dx * dx;
	    if(f1) (*n)->addForce(i, _viscosity * dx); 
      }
    }
  }

  double viscosity() const {return _viscosity;}

private:

  double _viscosity;
  blitz::Array<double, 1> _reference;
};


int main(int argc, char* argv[])
{
  std::cout << omp_get_max_threads() << " OpenMP threads." << std::endl;

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

  //double gamma_inp=0.0;
  double indent_inp=0.0;
  double viscosity_inp=0.0;
  double minStep_inp=0.01;
  double step_inp=0.0;
  bool prestressFlag = true;
  bool boundaryConditionFlag = false;
  bool roundAfmFlag = false;
  bool badCommandLine=false;
  for(char option; (option=getopt(argc,argv,"b:i:m:p:r:s:v:")) != EOF; ) 
    switch (option) {
    case 'b' :
      if( std::string(optarg) == std::string("yes") ) {
	boundaryConditionFlag = true;
      } else if( std::string(optarg) != std::string("no") ) {
	std::cout << "Yes or no for boundary conditions.  Exiting."
		  << std::endl;
	return 0;
      }	
      std::cout << "boundaryConditionFlag: " << boundaryConditionFlag
		<< std::endl;
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
      if( std::string(optarg) == std::string("no") ) 
	prestressFlag = false;
      std::cout << "prestressFlag: " << prestressFlag << std::endl;
      break;
    case 'r':
      roundAfmFlag = true;
      std::cout << "Round AFM." << std::endl;
      break;
    case 's':
      step_inp = std::atof(optarg);
      std::cout << "step size: " << step_inp << std::endl;
      break;
    case 'v':
      viscosity_inp = std::atof(optarg);
      std::cout << "viscosity: " << viscosity_inp << std::endl;
      break;
    default :
      badCommandLine=true;
      break;
    }

  if( badCommandLine || argc < 2 ) {
      cout << "Usage: indent modelName [-b boundaryConditionFlag -g gamma -p prestressFlag -n nsteps]." << endl;
      return(0);
  }
//   if(gamma_inp <=0.0) {
//     std::cout << "gamma should be positive." << std::endl;
//     return 0;
//   }

  string inputFileName = modelName;
  // create input stream
  ifstream ifs( inputFileName.c_str());
  if (!ifs) {
    cout << "Cannot open input file: " << inputFileName << endl;
    exit(0);
  }

  //
  // create vector of nodes
  double Rcapsid = 1.0;
  int dof=0;
  std::vector< NodeBase* > nodes;
  std::vector< XCNode<3>* > xcNodes;
  char key;  ifs>>key;
  double Ravg = 0;

  // find points header
  std::string token;
  ifs >> token; 
  while( token != "POINTS" ) ifs >> token;
  int npts=0;
  ifs >> npts; 
  xcNodes.reserve(npts);
  ifs >> token;// skip number type

  double eta = 0.0; //reaction coordinate
  srand(time(0));  
  // read in points
  for(int i=0; i<npts; i++) {
    int id=i;
    XCNode<3>::Point x;
    ifs >> x(0) >> x(1) >> x(2);
    Ravg += tvmet::norm2(x);
    NodeBase::DofIndexMap idx(4);
    for(int j=0; j<4; j++) idx[j]=dof++;
    eta= (double)(rand())/RAND_MAX;
    XCNode<3>* n = new XCNode<3>(id,idx,x,eta);
    nodes.push_back( n );
    xcNodes.push_back( n );
  }
  Ravg /= nodes.size();
  cout << "Number of nodes: " <<nodes.size() << endl
       << "Ravg = " << Ravg << endl;

  // read in triangle connectivities
  while( token != "POLYGONS" ) ifs >> token;
  vector< tvmet::Vector<int,3> > connectivities;
  tvmet::Vector<int, 3> c;
  int ntri=0; ifs >> ntri;
  connectivities.reserve(ntri);
  cout << "Number of triangles: " <<ntri << endl;
  int ntmp=0;
  ifs >> ntmp;
  if(ntmp != 4*ntri) {
    cout << "Non-triangular polygons?" << endl;
    return 0;
  }
  for (int i = 0; i<ntri; i++){
    int tmp=0;
    ifs >> tmp;
    ifs >> tmp; c[0]=tmp;
    ifs >> tmp; c[1]=tmp;
    ifs >> tmp; c[2]=tmp;
//     if(verbose) 
//       cout << setw(10) << c[0] << setw(10) << c[1] << setw(10) << c[2] << endl;
    connectivities.push_back(c);
  }
//   if(verbose) cout << connectivities.size() << endl;
  ifs.close();


  //double Y = sqrt(gamma_inp);
  double Y = 1.0;
  double KC = 1.0;
  double nu = 1.0/3.0;
  double KG = 0.0;
  double C0 = 0.0;
  double rho = 6.0;
  double gamma = 100.0;
  double DeltaEa = 0.01;
  double epsilon = 0.01;
  double DeltaEcPrime = 0.01;
  double Gamma = 1.0;
  double bdviscosity = 0.0;

  ifstream bdinp("constants.inp");
  bdinp >> rho >> gamma  >> nu >> epsilon >> Gamma >> DeltaEcPrime >> DeltaEa;
  std::cout << "rho ="     << rho << std::endl
            << "gamma ="   << gamma  <<  std::endl
            << "nu ="      << nu <<  std::endl
	    << "epsilon =" << epsilon    << std::endl
	    << "Gamma ="     << Gamma        << std::endl
	    << "DeltaEc' ="<< DeltaEcPrime   << std::endl
	    << "DeltaEa ="  << DeltaEa     << std::endl;

  double bb = -2.0-2.0*DeltaEcPrime;
  double cc = 1.0+3.0*DeltaEcPrime;
  double etaPrime = 0.5+1.5*DeltaEcPrime;
  double psiZero = DeltaEa/etaPrime/etaPrime/(etaPrime*etaPrime + bb*etaPrime + cc);
  double DeltaEc = DeltaEcPrime * psiZero;

  std::cout << "DeltaEc ="<< DeltaEc <<std::endl
	    << "psiZero ="<< psiZero <<std::endl;
  Y = gamma/rho/rho;
  std::cout << "Y ="<<Y<<std::endl;

  Rcapsid = rho;
  // rescale size 
  if( prestressFlag ) {
    for(int i=0; i<xcNodes.size(); i++) {
      XCNode<3>::Point x;
      x = xcNodes[i]->point();
      x *= Rcapsid/Ravg;
      xcNodes[i]->setPoint(x);
      xcNodes[i]->setPosition(x);
    }
  } else {
    // make capsid spherical 
    for(int i=0; i<xcNodes.size(); i++) {
      XCNode<3>::Point x;
      x = xcNodes[i]->point();
      double R = norm2(x);
      x *= Rcapsid/R;
      xcNodes[i]->setPoint(x);
      xcNodes[i]->setPosition(x);
    }
  }

  typedef GLElastic MaterialType;
  MaterialType protein(KC, KG, C0, Y, nu,  DeltaEa, epsilon, Gamma, DeltaEcPrime);
  std::cout << "GL Material has been created." << std::endl;
	
  // create Body
  int nBoundaries=0;
  double pressure=0.0;
  double tension=0.0;
  int quadOrder=2;
  double penaltyVolume=0.0;
  double penaltyArea=0.0;        
  GLBody<MaterialType> bd(protein, connectivities, nodes, quadOrder, nBoundaries, 
			  pressure, tension, penaltyVolume, penaltyArea, bdviscosity); 

  cout << "Created a body." << endl
       << "  # nodes: " << bd.nodes().size()<<std::endl
       << "  # shell elements: " << bd.elements().size() << std::endl;

  
  bd.setOutput(Body::paraview);

//   if( !prestressFlag ) {
//     bd.resetReference();
//     for(Body::ElementIterator e=bd.elements().begin(); e!=bd.elements().end(); e++) {
//       GLElement<MaterialType> * lse = (GLElement<MaterialType>*)(*e);
//       for(GLElement<MaterialType>::QuadPointIterator p=lse->quadraturePoints().begin(); p!=lse->quadraturePoints().end(); p++) {
// 	p->material.setSpontaneousCurvature(2.0*( p->material.meanCurvature() ) );
//       }
//     }
//   }

  bd.compute(true,true,false);

  // create Model
  Model::BodyContainer bdc;
  bdc.push_back(&bd);

  Model model(bdc);

//   model.checkConsistency(true,false);
//   model.checkRank(model.dof()-6,true);
//   return 0;

  int m=5;
  double factr=0.1;//1.0e+7;
  double pgtol=1.0e-5;
  int iprint = -1;
  double pentol=1.0e-4;
//   ifstream lbfgsbinp("lbfgsb.inp");
//   lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
//   if(verbose) 
    std::cout << "Input iprint: " << iprint << std::endl
	      << "Input factr: " << factr << std::endl
	      << "Input pgtol: " << pgtol << std::endl
	      << "Input m: " << m << std::endl
	      << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

  std::cout << "Relaxing shape for gamma = " << gamma << std::endl
	    << "Energy = " << solver.function() << std::endl;
  
  model.print("BeforeBFGS");
  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);

  std::cout << "Before BFGS:"<< std::endl
	    << "Volume = "<< bd.volume() << std::endl
	    << "Area = "<< bd.area() << std::endl
	    << "Reduced Volume = " << vred << std::endl
	    << "Energy = " << solver.function() << std::endl
	    << "strainEnergy = " << bd.totalStrainEnergy() << std::endl;
  // relax initial shape
  solver.solve(&model);

  string fname = modelName;
  fname += ".relaxed1000";
  model.print(fname);
  v = bd.volume();
  a = bd.area(); 
  vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);

    std::cout  << "Volume = "<< bd.volume() << std::endl
	       << "Area = "<< bd.area() << std::endl
	       << "Reduced Volume = " << vred << std::endl
	       << "Energy = " << solver.function() << std::endl
	       << "strainEnergy = " << bd.totalStrainEnergy() << std::endl;

  tvmet::Vector<double,3> Xavg(0.0);
  for ( int i = 0; i<xcNodes.size(); i++){
    Xavg += xcNodes[i]->point();
  }
  Xavg /= nodes.size();
  Ravg = 0.0;
  for ( int i = 0; i<xcNodes.size(); i++){
    Ravg += tvmet::norm2( xcNodes[i]->point() - Xavg );
  }
  Ravg /= nodes.size();
  
  double dRavg2 = 0.0;
  for ( int i = 0; i<xcNodes.size(); i++){
    double dR = (tvmet::norm2( xcNodes[i]->point() - Xavg) - Ravg); 
    dRavg2 += dR*dR;
  }
  dRavg2 /= nodes.size();
  
  double gammaCalc = Y*Ravg*Ravg/KC;
  double aspherity = dRavg2/(Ravg*Ravg);
  
  std::cout << "gamma = " << gammaCalc << endl
	    << "aspherity = " << aspherity << endl;

  std::cout << "Compressing capsid." << std::endl;

  // find top and bottom of capsid
  double Zmin=std::numeric_limits<double>::max();
  double Zmax=-std::numeric_limits<double>::max();
  double Zavg = 0.0;
  for (int a=0; a<xcNodes.size(); a++){
    double Z = xcNodes[a]->getPoint(2);
    Zmin = std::min(Zmin, Z);
    Zmax = std::max(Zmax, Z);
    Zavg += Z;
  }
  Zavg /= xcNodes.size();

  // create indentor and plate
  double afmR = rho;
  tvmet::Vector<double,3> xc; xc = 0.0, 0.0, Zmax+afmR; 

  RigidHemisphereContactSPT * afm 
    = new RigidHemisphereContactSPT( xcNodes, afmR, xc );
  if(roundAfmFlag) bd.pushBackContact( afm );

  // set up bounds for solver
  blitz::Array<int,1> nbd(4*nodes.size());
  blitz::Array<double,1> lo(4*nodes.size());
  blitz::Array<double,1> hi(4*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;
  for (int a=0; a<xcNodes.size(); a++) {  
    if( roundAfmFlag )  nbd(4*a+2) = 1;
    else 		nbd(4*a+2) = 2;
    hi(4*a+2) = Zmax;
    lo(4*a+2) = Zmin;
  }
  if(boundaryConditionFlag) {
    // constrain apices
    for (int a=0; a<xcNodes.size(); a++) {
      double x = xcNodes[a]->getPoint(0);
      double y = xcNodes[a]->getPoint(1);
      double R = sqrt(x*x + y*y);
      //std:: cout << R << std::endl;
      if(R < 1.0e-5) {
	std::cout << "Constraining node " << a << " at R=" << R 
		  << "." << std::endl;
	for(int i=0; i<2; i++) {
	  nbd(4*a+i) = 2;
	  hi(4*a+i) = 0.0;
	  lo(4*a+i) = 0.0;
	}
      }
    }
  }

  // add some viscosity for regularization
  double viscosity=viscosity_inp;
  ViscousRegularizer vr(bd.nodes(), viscosity);
  bd.pushBack( &vr ); 

  double vrTol = 1.0e-6;

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
  double vrEnergy = vr.energy();
  double bdEnergy = bd.energy();

  string fzName = modelName + ".fz";
  ofstream FvsZ(fzName.c_str());

  double F_prev = -1.0;
  double Z_drop = Zbegin;
  blitz::Array<double,1> x_prev(model.dof());
  model.getField(solver);
  for(int i=0; i<model.dof(); i++ ) x_prev(i) = solver.field(i);

  for(double Z = Zbegin; Z<Zend; Z+=dZ) {
    for (int a=0; a<xcNodes.size(); a++) {  
      lo(4*a+2) = Z;
    }
    solver.setBounds(nbd,lo,hi);

    int viter=0;
    do {
      if(verbose) std::cout << "viter = " << viter << std::endl;
      viter++;
      vr.step();
      solver.solve( &model );
      vrEnergy = vr.energy();
      bdEnergy = bd.energy();
#ifdef WITH_MPI
      double temp=vrEnergy;
      MPI_Allreduce(&temp, &vrEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
      temp=bdEnergy;
      MPI_Allreduce(&temp, &bdEnergy, 1, MPI_DOUBLE, 
		    MPI_SUM, MPI_COMM_WORLD);
#endif      
      if(verbose)
	std::cout << "v: " << vrEnergy << std::endl
		  << "e: " << bdEnergy << std::endl;
    } while(vrEnergy > vrTol);
    double height = Zmax-Z;
    // add up forces on top and bottom
    double FZg = 0.0;
    double FZa = 0.0;
    double Ztop = Zmax;
    double Zbot = Zmax - height;
    for (int a=0; a<xcNodes.size(); a++) {
      double F = xcNodes[a]->getForce(2);
      if( std::abs( xcNodes[a]->getPoint(2) - Ztop ) < 1.0e-4 ) {
	// top
	FZa += xcNodes[a]->getForce(2);
      } else if( std::abs( xcNodes[a]->getPoint(2) - Zbot ) < 1.0e-4 ) {
	// bottom
	FZg += xcNodes[a]->getForce(2);
      }
    }
    if(roundAfmFlag) FZa = afm->FZ();

#ifdef WITH_MPI
    temp=FZg;
    MPI_Reduce(&temp, &FZg, 1, MPI_DOUBLE, 
	       MPI_SUM, 0, MPI_COMM_WORLD);
    temp=FZa;
    MPI_Reduce(&temp, &FZa, 1, MPI_DOUBLE, 
	       MPI_SUM, 0, MPI_COMM_WORLD);
#endif      
    if( FZg < F_prev && dZ > minStep_inp ) {
      // drop in force.  Back up and try a smaller increment
      Z_drop = Z;
      std::cout << "Force drop after increment of "<< dZ << std::endl;
      for(int i=0; i<model.dof(); i++ ) solver.field(i) = x_prev(i);
      model.putField(solver);
      Z -= dZ;
      dZ /= 2.0;
      std::cout << "Trying again with increment of " << dZ << std::endl;
      continue;
    } else if ( Z > Z_drop && dZ < step_inp ) {
      dZ *= 2.0;
      std::cout << "Increasing increment to " << dZ << std::endl;
    }
    F_prev = FZg;

    char name[20]; 
    sprintf(name,"1000disp%f",originalHeight-height);
    model.print(name);
    if(verbose) {
      std::cout << "Energy = " << solver.function() << std::endl
		<< "  height = " << height
		<< "  change in height = " << originalHeight-height
		<< "  FZg = " << FZg
		<< "  FZa = " << FZa << std::endl;
    }
    FvsZ << std::setw( 24 ) << std::setprecision(16) << originalHeight-height
	 << std::setw( 24 ) << std::setprecision(16) << FZg
	 << std::setw( 24 ) << std::setprecision(16) << FZa
	 << std::setw( 24 ) << std::setprecision(16) << solver.function()
	 << std::endl;
  }
  FvsZ.close();

  std::cout << "All done.  Bye now." << std::endl;

  return 0;
}

