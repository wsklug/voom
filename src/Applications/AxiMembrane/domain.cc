#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "FVK.h"
#include "C1AxiShellBody.h"
#include "ConcentratedForce.h"
#include "ViscousRegularizer.h"
#include "Spring.h"
#include "Model.h"
#include "Lbfgsb.h"

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{
  typedef EvansElastic Material_type;

  typedef C1AxiShellBody<Material_type> Body_type;

  // make a straight line of nodes and elements, from r=0 to r=R
  int N0 = 15;
  double R0 = 0.5;
  int nDomainNodes = 2*(N0+1);
  int multiple=6;
  double R = multiple*R0;
  int N = multiple*N0;
  double dr = R/N;

  int nAllNodes = 2*(N+1);
  std::vector< NodeBase* > allNodes;
  allNodes.reserve(nAllNodes);

  int dof=0;
  for(int i=0; i<=N; i++) {
    double r=i*dr;
    DeformationNode<2>::Point X;
    DeformationNode<2>::Point x;
    X = r, 0.0;
    x = X;
    NodeBase::DofIndexMap idx(2);
    idx[0]=dof++; idx[1]=dof++;
    allNodes.push_back(new DeformationNode<2>(2*i,idx,X,x));    
    
    DeformationNode<2>::Point v;
    v = dr, 0.0;
    idx[0]=dof++; idx[1]=dof++;
    allNodes.push_back(new DeformationNode<2>(2*i+1,idx,v));        
  }

  Body_type::ConnectivityContainer domainConnect;
  for(int i=0; i<N0; i++) {
    Body_type::ElementConnectivity c;
    c = 2*i, 2*i+1, 2*i+2, 2*i+3;
    domainConnect.push_back( c );
  }

  Body_type::ConnectivityContainer membraneConnect;
  for(int i=N0; i<N; i++) {
    Body_type::ElementConnectivity c;
    c = 2*i, 2*i+1, 2*i+2, 2*i+3;
    membraneConnect.push_back( c );
  }

  cout << allNodes.size() << endl;
  for(int i=0; i<nAllNodes; i++) {
    cout << setw(12) << allNodes[i]->id()
	 << setw(20) << allNodes[i]->getPoint(0)
	 << setw(20) << allNodes[i]->getPoint(1)
	 << endl;
  }	

  cout << "Domain connectivities " << domainConnect.size() << endl;
  for(int i=0; i<domainConnect.size(); i++) {
    cout << setw(12) << domainConnect[i][0] 
	 << setw(12) << domainConnect[i][1] 
	 << setw(12) << domainConnect[i][2]
	 << setw(12) << domainConnect[i][3]
	 << endl;
  }

  cout << "Membrane connectivities " << membraneConnect.size() << endl;
  for(int i=0; i<membraneConnect.size(); i++) {
    cout << setw(12) << membraneConnect[i][0] 
	 << setw(12) << membraneConnect[i][1] 
	 << setw(12) << membraneConnect[i][2]
	 << setw(12) << membraneConnect[i][3]
	 << endl;
  }

  // no pressure
  double pressure=0.0;

  // do displacement control, set applied tension to zero
  double tension=1.0;

  // no global constraints
  double volumeModulus=0.0;
  double areaModulus=1.0e3;

  unsigned quadOrder = 5;
  
  double KC=1.0;
  double KG=0.0;
  double C0=0.0;
  double KA=1.0e+3;
  double G = 0.0;
  Material_type bilayer(KC,KG,C0,G,KA);
  std::cout << "EvansElastic Material has been created." << std::endl;

  C1AxiShellBody<Material_type> membraneBody(bilayer, 
					     membraneConnect, 
					     allNodes, 
					     quadOrder, 
					     pressure, 
					     tension,
					     0.0,
					     volumeModulus, 
					     areaModulus,
					     0.0,
					     noConstraint,
					     augmented //multiplier
					     );
    
  C1AxiShellBody<Material_type> domainBody(bilayer, 
					   domainConnect, 
					   allNodes, 
					   quadOrder, 
					   pressure, 
					   tension,
					   0.0,
					   volumeModulus, 
					   areaModulus,
					   0.0,
					   noConstraint,
					   augmented // augmented
					   );

  // add line tension to interface node
  double gamma=0.0;
  Vector2D LT;
  LT = -2.0*M_PI*gamma, 0.0;
  DeformationNode<2> * interface 
    = static_cast< DeformationNode<2>* >(allNodes[2*N0]);
  ConcentratedForce<2> lineTension(interface, LT);
  domainBody.addElement( &lineTension );


  // add some viscous regularization for buckling
  double initialViscosity = 1.0e0;
  double minViscosity = 1.0e-3;
  ViscousRegularizer reg(allNodes, initialViscosity);

  // try without viscosity
  // domainBody.addElement( &reg );

  // add some springs to regularize element lengths
  double springConstant=1.0e4;
  std::vector<Spring<2>*> springs;
  for(int I=0; I<nDomainNodes-2; I+=2) {
    DeformationNode<2> * ndA = dynamic_cast< DeformationNode<2>* >(allNodes[I]);
    DeformationNode<2> * ndB = dynamic_cast< DeformationNode<2>* >(allNodes[I+2]);
    assert( ndA != 0 );
    assert( ndB != 0 );  
    Spring<2> * sp = new Spring<2>(ndA, ndB, springConstant);
//     springs.push_back( sp );
//     domainBody.addElement( sp );
  }
  
  springConstant = 1.0e2;
  for(int I=nDomainNodes; I<nAllNodes-2; I+=2) {
    DeformationNode<2> * ndA = dynamic_cast< DeformationNode<2>* >(allNodes[I]);
    DeformationNode<2> * ndB = dynamic_cast< DeformationNode<2>* >(allNodes[I+2]);
    assert( ndA != 0 );
    assert( ndB != 0 );  
    Spring<2> * sp = new Spring<2>(ndA, ndB, springConstant);
//     springs.push_back( sp );
//     domainBody.addElement( sp );
  }

  Model::BodyContainer bodies;
  bodies.push_back( &domainBody );
  bodies.push_back( &membraneBody);

  Model model( bodies, allNodes );

//   model.checkConsistency(true,false);
//   model.print("check");
//   model.checkRank(2*nMembraneNodes-1, true);
//   return 0;

  int m=5;
  double factr=1.0e+7;
  double pgtol=1.0e-5;
  int iprint = 0;
  double pentol=1.0e-4;
  ifstream lbfgsbinp("lbfgsb.inp");
  lbfgsbinp >> iprint >> factr >> pgtol >> m >> pentol;
  std::cout << "Input iprint: " << iprint << std::endl
	    << "Input factr: " << factr << std::endl
	    << "Input pgtol: " << pgtol << std::endl
	    << "Input m: " << m << std::endl
	    << "Input pentol: " << pentol << std::endl;
  Lbfgsb solver(model.dof(), m, factr, pgtol, iprint );//(true);

  // set up bounds for solver
  blitz::Array<int,1> nbd(2*allNodes.size());
  blitz::Array<double,1> lo(2*allNodes.size());
  blitz::Array<double,1> hi(2*allNodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;

  // keep all position nodes at positive r
  for(int I=0; I<nAllNodes; I+=2) {
    nbd(2*I+0) = 1;
    lo(2*I+0) = 0.0;
  }


//   domainBody.resetReference();
//   membraneBody.resetReference();

  // first node, fix r-position
  int I=0;
  nbd(2*I+0) = 2;
  hi(2*I+0) = allNodes[I]->getPoint(0);
  lo(2*I+0) = allNodes[I]->getPoint(0);
  
  // first node, fix tangent
  I=1;
  for(int i=0; i<2; i++) {
    nbd(2*I+i) = 2;
    hi(2*I+i) = allNodes[I]->getPoint(i);
    lo(2*I+i) = allNodes[I]->getPoint(i);
  }

  // interface node, fix z-position
  I=2*N0;
  nbd(2*I+1) = 2;
  hi(2*I+1) = allNodes[I]->getPoint(1);
  lo(2*I+1) = allNodes[I]->getPoint(1);

  // last node, fix r-position
  I=nAllNodes-2;
  nbd(2*I+0) = 2; // 1 for lower bound only
  hi(2*I+0) = allNodes[I]->getPoint(0);
  lo(2*I+0) = allNodes[I]->getPoint(0);

  solver.setBounds(nbd,lo,hi);

  std::cout << "Initial state:"<<std::endl
	    << std::endl
	    << "CONSTRAINTS:" << std::endl
	    << "     Domain Area = "<< domainBody.area() << std::endl
	    << "    (prescribed) = "<< domainBody.prescribedArea() << std::endl
	    << "  Domain tension = " << domainBody.tension() << std::endl
	    << "   Membrane Area = "<< membraneBody.area() << std::endl
	    << "Membrane tension = " << membraneBody.tension() << std::endl
	    << std::endl
	    << "ENERGY:" << std::endl
	    << "     Total = " << solver.function() << std::endl
	    << "   Bending = " << domainBody.bendingEnergy() << std::endl
	    << "Stretching = " << domainBody.stretchingEnergy() << std::endl
	    << " Interface = " << lineTension.energy() << std::endl
	    << "   Viscous = " << reg.energy() << std::endl
	    << std::endl
	    << "DEFORMED POSITIONS:" << std::endl
	    << "first node: " 
	    << std::setw(20) << allNodes[0]->getPoint(0)
	    << std::setw(20) << allNodes[0]->getPoint(1)
	    << std::endl
	    << " interface: " 
	    << std::setw(20) << interface->getPoint(0)
	    << std::setw(20) << interface->getPoint(1)
	    << std::endl
	    << " last node: " 
	    << std::setw(20) << allNodes[nAllNodes-2]->getPoint(0)
	    << std::setw(20) << allNodes[nAllNodes-2]->getPoint(1)
	    << std::endl
	    << "============================================================="
	    << std::endl;

  model.computeAndAssemble(solver, true, true, false);
  model.print("initial"); 
  double a = domainBody.area(); 
  double A = domainBody.prescribedArea();

  char fname[20];
 
  
//   for(gamma = 15.0; gamma <= 20.0; gamma += 1.0) {
  gamma = 13.14;
  lineTension.setForce(0,-2.0*M_PI*gamma);

  domainBody.setFixedTension( tension - gamma/R0 );
  
  double delta_r = -2.0e-2*dr;
  int n_steps = 100;
  double r_out =  allNodes[I]->getPoint(0);

//   r_out -= (n_steps-3)*delta_r;
  hi(2*I+0) = 
  lo(2*I+0) = r_out;
  solver.setBounds(nbd,lo,hi);
  

//   r_out -= delta_r;
//   hi(2*I+0) = lo(2*I+0) = r_out;


  //perturb all position nodes in domain up
  for(int i=0; i<N0; i+=2) {
    double r=allNodes[i]->getPoint(0);
    double amplitude=1.0e-1;
    double z=amplitude*cos(0.5*M_PI*r/R0);
//     double dz=-amplitude*dr*(0.5*M_PI/R0)*sin(0.5*M_PI*r/R0);
    allNodes[i]->setPoint(1, z);
//     allNodes[i+1]->setPoint(1, dz);
  }
  // point each director toward next node
  for(int i=0; i<N0; i+=2) {
    double dz=allNodes[i+2]->getPoint(1) - allNodes[i]->getPoint(1);
    allNodes[i+1]->setPoint(1,dz);
  }  

  ofstream ofs("domain.dat");

  for(int step=0; step<n_steps; step++) {

    std::cout << std::endl
	      << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	      << std::endl
	      << " step = " << step
	      << " gamma = " << gamma 
	      << " r_out = " << r_out 
	      << std::endl
	      << std::endl;

//     //perturb nodes
//     for(int i=0; i<nAllNodes; i+=2) {
//       double r=allNodes[i]->getPoint(0);
//       double amplitude=1.0e-1;
//       double z=amplitude*cos(0.5*M_PI*r/R0);
//       double dz=-amplitude*dr*(0.5*M_PI/R0)*sin(0.5*M_PI*r/R0);
//       allNodes[i]->addPoint(1, z);
//       allNodes[i+1]->addPoint(1, dz);
//     }

    areaModulus=1.0e3;
    double viscosity=initialViscosity;     
    double targetVelocity = 1.0e-2;
    for(int iter=0; iter<100; iter++) {
      std::cout << "ITERATION "<< iter << ":"
		<< std::endl 
		<< std::endl;
      reg.step();
//       membraneBody.resetReference();
//       domainBody.resetReference();
      solver.solve(&model);
      sprintf(fname,"iter-%02d",iter);
      //model.print(fname);

      domainBody.updateFixedTension();
      membraneBody.updateFixedTension();
      domainBody.updatePenaltyArea(areaModulus*=3);
      membraneBody.updatePenaltyArea(areaModulus*=3);
    
      a = domainBody.area(); 
      A = domainBody.prescribedArea(); 
      
      double springEnergy = 0.0;
      for(int I=0; I<springs.size(); I++) {
	springEnergy += springs[I]->energy();
	springs[I]->resetLength();
      }

      viscosity *= reg.velocity()/targetVelocity;
      viscosity = std::max(viscosity,minViscosity);
      reg.setViscosity(viscosity);

      std::cout << "CONSTRAINTS:" << std::endl
		<< std::setprecision(16) 
		<< "     Domain Area = "<< domainBody.area() << std::endl
		<< "    (prescribed) = "<< domainBody.prescribedArea() << std::endl
		<< "  Domain tension = " << domainBody.tension() << std::endl
		<< "   Membrane Area = "<< membraneBody.area() << std::endl
		<< "    (prescribed) = "<< membraneBody.prescribedArea() << std::endl
		<< "Membrane tension = " << membraneBody.tension() << std::endl
		<< std::endl
		<< std::setprecision(6) 
		<< "ENERGY:" << std::endl
		<< "     Total = " << solver.function() << std::endl
		<< "   Bending = " << domainBody.bendingEnergy() << std::endl
		<< "Stretching = " << domainBody.stretchingEnergy() << std::endl
		<< " Interface = " << lineTension.energy() << std::endl
		<< "   Viscous = " << reg.energy() << std::endl
		<< "   Springs = " << springEnergy << std::endl
		<< std::endl
		<< "DEFORMED POSITIONS:" << std::endl
		<< "first node: "
		<< std::setw(20) << allNodes[0]->getPoint(0)
		<< std::setw(20) << allNodes[0]->getPoint(1)
		<< std::endl
		<< " interface: " 
		<< std::setw(20) << interface->getPoint(0)
		<< std::setw(20) << interface->getPoint(1)
		<< std::endl
		<< " last node: " 
		<< std::setw(20) << allNodes[nAllNodes-2]->getPoint(0)
		<< std::setw(20) << allNodes[nAllNodes-2]->getPoint(1)
		<< std::endl
		<< std::endl
		<< "VISCOSITY: " << std::endl
		<< " velocity = " << reg.velocity() << std::endl
		<< " updated viscosity = " << viscosity << std::endl
		<< std::endl
		<< "============================================================="
		<< std::endl;

      if( reg.energy()+springEnergy < 1.0e-8*domainBody.energy() 
	  && abs(domainBody.area() - domainBody.prescribedArea()) < 1.0e-6 )
//       if( reg.energy()+springEnergy < 1.0e-8*domainBody.energy()  ) 
	break;

    }

    //sprintf(fname,"gamma-%f",gamma);
    sprintf(fname,"step-%03d",step);
    model.print(fname);

    ofs << std::setprecision(16) 
	<< std::setw(24) << r_out 
	<< std::setw(24) << domainBody.tension() 
	<< std::setw(24) << membraneBody.tension() 
	<< std::setw(24) << allNodes[0]->getPoint(0)
	<< std::setw(24) << allNodes[0]->getPoint(1)
	<< std::setw(24) << interface->getPoint(0)
	<< std::setw(24) << interface->getPoint(1)
	<< std::setw(24) << allNodes[nAllNodes-2]->getPoint(0)
	<< std::setw(24) << allNodes[nAllNodes-2]->getPoint(1)
	<< std::endl; 
      


    r_out += delta_r;

    hi(2*I+0) = 
    lo(2*I+0) = r_out;
    solver.setBounds(nbd,lo,hi);
    

    // if(membraneBody.tension() <= 1.0) break;

  }
  
  ofs.close();
  return 0;
}

