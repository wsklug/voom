#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "FVK.h"
#include "StVenantShell.h"
#include "C1AxiShellBody.h"
#include "Model.h"
#include "Lbfgsb.h"

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{

  double R = 500.0; // nm

  double KC=20.0; // kBT
  double KG=-0.5*KC;
  double C0=0.0;

  typedef SCElastic MembraneMaterial;
  MembraneMaterial bilayer(KC,KG,C0);
  std::cout << "SCElastic Material has been created." << std::endl;

  
  double E = 50.0/4.1; // kBT/nm^3 
  double h = 5.0;      // nm
  double Y = E*h;
  double nu=0.5;
  KC = E*pow(h,3)/(12.0*(1-sqr(nu)));
  KG = -2.0*(1.0-nu)*KC;
  C0 = -1.0/R;
  typedef StVenantShell WallMaterial;
  WallMaterial peptidoglycan(E,nu,h);
  std::cout << "FVK Material has been created." << std::endl;

  // double KA=1.0e2;
  // double G = 0.0;
  // typedef EvansElastic Material_type;
  // Material_type bilayer(KC,KG,C0,G,KA);
  // std::cout << "EvansElastic Material has been created." << std::endl;

  typedef C1AxiShellBody<MembraneMaterial> MembraneBody;
  typedef C1AxiShellBody<WallMaterial> WallBody;

  MembraneBody::ConnectivityContainer connectivities;

  // make a half-circle of nodes and elements, from -Pi/2 to Pi/2.
  double L= R;
  int NL = 100;
  int NR = 130;
  int N = NL + NR;
  int nNodes = 2*(N+1);
  int nElements = N;

  std::vector< NodeBase* > nodes;
  nodes.reserve(nNodes);

  double dtheta = M_PI/NR;
  double dL=L/NL;

  int dof=0;
  //
  // points for cylindrical section
  //
  for(  int i=0; i<NL; i++) {
    double z = i*dL;
    DeformationNode<2>::Point X;
    DeformationNode<2>::Point x;
    X = R, z;
    x = X(0), X(1); // stretch a little in the z-dir?
    NodeBase::DofIndexMap idx(2);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*i,idx,X,x));    

    DeformationNode<2>::Point v;
    v = 0.0, dL;
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*i+1,idx,v));        

    MembraneBody::ElementConnectivity c;
    c = 2*i, 2*i+1, 2*i+2, 2*i+3;
    connectivities.push_back( c );
  }
  //
  // points for end cap
  //
  for(  int i=0; i<NR; i++) {
    double theta=i*0.5*M_PI/NR;
    DeformationNode<2>::Point X;
    DeformationNode<2>::Point x;
    X = R*cos(theta), L+R*sin(theta);
    x = X(0), X(1); // stretch a little in the z-dir?
    NodeBase::DofIndexMap idx(2);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*(i+NL),idx,X,x));    

    DeformationNode<2>::Point v;
    v = -R*dtheta*sin(theta), R*dtheta*cos(theta);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*(i+NL)+1,idx,v));        

    MembraneBody::ElementConnectivity c;
    c = 2*(i+NL), 2*(i+NL)+1, 2*(i+NL)+2, 2*(i+NL)+3;
    connectivities.push_back( c );
  }
  {
    DeformationNode<2>::Point x;
    x = 0.0, L+R;
    NodeBase::DofIndexMap idx(2);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*N,idx,x));    
    
    DeformationNode<2>::Point v;
    v = -R*dtheta, 0.0;
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*N+1,idx,v));        
  }

  cout << nodes.size() << endl;
  for(int i=0; i<nNodes; i++) {
    cout << setw(12) << nodes[i]->id()
	 << setw(20) << nodes[i]->getPoint(0)
	 << setw(20) << nodes[i]->getPoint(1)
	 << endl;
  }	

  cout << connectivities.size() << endl;
  for(int i=0; i<connectivities.size(); i++) {
    cout << setw(12) << connectivities[i][0] 
	 << setw(12) << connectivities[i][1] 
	 << setw(12) << connectivities[i][2]
	 << setw(12) << connectivities[i][3]
	 << endl;
  }
  
	
  //
  // create Body
  double atm = 0.025;
  double pressure=10.0*atm; // 1atm = 0.025 kBT/nm^2
  double tension=0;
  double volumeModulus=0.0;//1.0e3;
  double areaModulus=60.0; // kBT/nm^2

  unsigned quadOrder = 5;
  
  //
  // Membrane body
  //
  MembraneBody membrane(bilayer, 
			connectivities, 
			nodes, 
			quadOrder, 
			pressure, 
			tension,
			0.0,
			volumeModulus, 
			areaModulus,
			0.0,
			noConstraint,
			penalty
			);
  
  cout << "Created a body." << endl;

  membrane.compute(true,true, false);
  cout << "volume of current body = " << membrane.volume() << endl;
  cout << "energy of current body = " << membrane.energy() << endl;

  Model::BodyContainer bodies;
  bodies.push_back( &membrane );

  //
  // Wall body
  //
  WallBody wall(peptidoglycan, 
		connectivities, 
		nodes, 
		quadOrder, 
		0.0, 
		0.0,
		0.0,
		volumeModulus, 
		areaModulus,
		0.0,
		noConstraint,
		noConstraint
		);
  
  cout << "Created a body." << endl;

  wall.compute(true,true, false);
  cout << "volume of current body = " << wall.volume() << endl;
  cout << "energy of current body = " << wall.energy() << endl;

  bodies.push_back( &wall );

  Model model( bodies, nodes );

//   for(int I=0; I<nNodes; I++) {
//     double x=nodes[I]->getPoint(0);
//     x *= 0.75;
//     nodes[I]->setPoint(0,x);
//   }

//   bd.compute(true,true, false);
//   cout << "volume of current body = " << bd.volume() << endl;
//   cout << "energy of current body = " << bd.energy() << endl;

//   model.checkConsistency(true,false);
//   model.print("check");
//   model.checkRank(2*nNodes-1, true);
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
  blitz::Array<int,1> nbd(2*nodes.size());
  blitz::Array<double,1> lo(2*nodes.size());
  blitz::Array<double,1> hi(2*nodes.size());
  nbd = 0;
  lo = 0.0;
  hi = 0.0;

  // keep all position nodes at positive r
  for(int I=0; I<nNodes; I+=2) {
    nbd(2*I+0) = 1;
    lo(2*I+0) = 0.0;
  }

  // first node
  int I=0;
  // // r-position
  // nbd(2*I+0) = 2;
  // hi(2*I+0) = R;
  // lo(2*I+0) = R;

  // z-position
  nbd(2*I+1) = 2;
  hi(2*I+1) = 0.0;
  lo(2*I+1) = 0.0;

  // second node: should be tangent to z; constratin r
  I=1;
  // r-tangent
  nbd(2*I+0) = 2;
  hi(2*I+0) = nodes[I]->getPoint(0);
  lo(2*I+0) = nodes[I]->getPoint(0);
  // z-tangent
  //nbd(2*I+1) = 2;
  //hi(2*I+1) = nodes[I]->getPoint(1);
  //lo(2*I+1) = nodes[I]->getPoint(1);

  // second to last node
  I=nNodes-2;
  // r-position
  nbd(2*I+0) = 2;
  hi(2*I+0) = 0.0;
  lo(2*I+0) = 0.0;

  // last node: should be tangent to r; constrain z
  I=nNodes-1;
  // r-tangent
  //nbd(2*I+0) = 2;
  //hi(2*I+0) = nodes[I]->getPoint(0);
  //lo(2*I+0) = nodes[I]->getPoint(0);
  // z-tangent
  nbd(2*I+1) = 2;
  hi(2*I+1) = nodes[I]->getPoint(1);
  lo(2*I+1) = nodes[I]->getPoint(1);

  solver.setBounds(nbd,lo,hi);

  model.print("initial"); 
  double v = membrane.volume();
  double a = membrane.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = membrane.prescribedVolume(); 
  double A = membrane.prescribedArea();
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);

  ofstream output("output.dat");

  output.precision(8);

  for(int step = 0; step<=20; step++) {

    double p = pressure*step/20.0;
    membrane.setFixedPressure(p);

    char fname[20];

    solver.solve(&model);
    sprintf(fname,"membrane-%02d",step);
    // model.print(fname);
    membrane.printParaview(fname);
    sprintf(fname,"wall-%02d",step);
    // model.print(fname);
    wall.printParaview(fname);
      
    v = 2.0*membrane.volume();
    a = 2.0*membrane.area(); 
    vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
    V = 2.0*membrane.prescribedVolume();
    A = 2.0*membrane.prescribedArea(); 
    Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);

    double pLY=0.5*R*membrane.fixedPressure(); // sqrt(0.25*a/M_PI)

    std::cout << "step "<< step << ":"<< std::endl
	      << "Volume = "<< membrane.volume() << std::endl
	      << "Initial Volume = "<< membrane.prescribedVolume() << std::endl
	      << "Area = "<< membrane.area() << std::endl
	      << "Initial Area = "<< membrane.prescribedArea() << std::endl
	      << "Reduced Volume = " << vred << std::endl
	      << "Initial Reduced Volume = " << Vred << std::endl
	      << "Volume change = " << v/V << std::endl
	      << "Area change   = " << a/A << std::endl
	      << "Energy = " << membrane.energy() << std::endl
	      << "Stretching Energy = " << membrane.stretchingEnergy() << std::endl
	      << "pressure = " << membrane.fixedPressure() << std::endl
	      << "tension = " << membrane.tension() << std::endl
	      << "Laplace tension = " << pLY << std::endl;

    output << std::setw(8)  << step
	   << std::setw(14) << membrane.fixedPressure()
	   << std::setw(14) << membrane.tension()
	   << std::setw(14) << a
	   << std::setw(14) << a/A
	   << std::setw(14) << v
	   << std::setw(14) << v/V 
	   << std::endl;
  }
  output.close();

  return 0;
}

