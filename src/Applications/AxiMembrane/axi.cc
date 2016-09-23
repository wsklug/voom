#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "Node.h"
#include "SCElastic.h"
#include "EvansElastic.h"
#include "FVK.h"
#include "C1AxiShellBody.h"
#include "Model.h"
#include "Lbfgsb.h"

using namespace tvmet;
using namespace std;
using namespace voom;

void ioSetting(int argc, char* argv[], ifstream&, string&);

int main(int argc, char* argv[])
{
  double KC=1.0;
  double KG=0.0;
  double C0=0.0;

//   typedef SCElastic Material_type;
//   Material_type bilayer(KC,KG,C0);
//   std::cout << "SCElastic Material has been created." << std::endl;

//   double Y=1.0e5;
//   double nu=0.3;
//   typedef FVK Material_type;
//   Material_type bilayer(KC,KG,C0,Y,nu);
//   std::cout << "FVK Material has been created." << std::endl;

  double KA=1.0e2;
  double G = 0.0;
  typedef EvansElastic Material_type;
  Material_type bilayer(KC,KG,C0,G,KA);
  std::cout << "EvansElastic Material has been created." << std::endl;

  typedef C1AxiShellBody<Material_type> Body_type;


  Body_type::ConnectivityContainer connectivities;

  // make a half-circle of nodes and elements, from -Pi/2 to Pi/2.
  int N = 50;
  int nNodes = 2*(N+1);
  int nElements = N;

  std::vector< NodeBase* > nodes;
  nodes.reserve(nNodes);

  double R = 1.0;
  double dtheta = M_PI/N;
  double h= 0.0;

  int dof=0;
  for(  int i=0; i<N; i++) {
    double theta=i*M_PI/N - M_PI/2;
    DeformationNode<2>::Point X;
    DeformationNode<2>::Point x;
    X = h + R*cos(theta), R*sin(theta);
    x = X(0), X(1); // stretch a little in the z-dir?
    NodeBase::DofIndexMap idx(2);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*i,idx,X,x));    

    DeformationNode<2>::Point v;
    v = -R*dtheta*sin(theta), R*dtheta*cos(theta);
    idx[0]=dof++; idx[1]=dof++;
    nodes.push_back(new DeformationNode<2>(2*i+1,idx,v));        

    Body_type::ElementConnectivity c;
    c = 2*i, 2*i+1, 2*i+2, 2*i+3;
    connectivities.push_back( c );
  }
  {
    DeformationNode<2>::Point x;
    x = h, R;
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
  double pressure=-10.0;
  double tension=-5.0;
  double volumeModulus=1.0e3;
  double areaModulus=1.0e3;

  unsigned quadOrder = 5;
  
  C1AxiShellBody<Material_type> bd(bilayer, 
				   connectivities, 
				   nodes, 
				   quadOrder, 
				   pressure, 
				   tension,
				   0.0,
				   volumeModulus, 
				   areaModulus,
				   0.0,
				   augmented,
				   augmented
				   );

  cout << "Created a body." << endl;

  bd.compute(true,true, false);
  cout << "volume of current body = " << bd.volume() << endl;
  cout << "energy of current body = " << bd.energy() << endl;

  Model::BodyContainer bodies;
  bodies.push_back( &bd );

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
  // r-position
  nbd(2*I+0) = 2;
  hi(2*I+0) = 0.0;
  lo(2*I+0) = 0.0;
  // z-position
  nbd(2*I+1) = 2;
  hi(2*I+1) = nodes[I]->getPoint(1);
  lo(2*I+1) = nodes[I]->getPoint(1);

  // second node
  I=0;
  // r-tangent
  nbd(2*I+0) = 2;
  hi(2*I+0) = nodes[I]->getPoint(0);
  lo(2*I+0) = nodes[I]->getPoint(0);
  // z-tangent
  nbd(2*I+1) = 2;
  hi(2*I+1) = nodes[I]->getPoint(1);
  lo(2*I+1) = nodes[I]->getPoint(1);

  // second to last node
  I=nNodes-2;
  // r-position
  nbd(2*I+0) = 2;
  hi(2*I+0) = 0.0;
  lo(2*I+0) = 0.0;

  // last node
  I=nNodes-1;
  // r-tangent
  nbd(2*I+0) = 2;
  hi(2*I+0) = nodes[I]->getPoint(0);
  lo(2*I+0) = nodes[I]->getPoint(0);
  // z-tangent
  nbd(2*I+1) = 2;
  hi(2*I+1) = nodes[I]->getPoint(1);
  lo(2*I+1) = nodes[I]->getPoint(1);

  // squash in z-dir
  for(int I=0; I<nNodes; I+=2) {
      nbd(2*I+1) = 2;
      hi(2*I+1) =  0.9*R;
      lo(2*I+1) = -0.9*R;
  }
  solver.setBounds(nbd,lo,hi);

  model.print("initial"); 
  double v = bd.volume();
  double a = bd.area(); 
  double vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
  double V = bd.prescribedVolume(); 
  double A = bd.prescribedArea();
  double Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);

  for(Vred=0.98; Vred>0.5; Vred-=0.02) {

    char fname[20];

    bd.reduceVolume(Vred);

    volumeModulus=1.0e3;
    areaModulus=1.0e3;
    for(int iter=0; iter<10; iter++) {
      bd.resetReference();
      solver.solve(&model);
      sprintf(fname,"iter-%02d",iter);
      model.print(fname);
      bd.updateFixedPressure();
      bd.updateFixedTension();
      bd.updatePenaltyVolume(volumeModulus*=3);
      bd.updatePenaltyArea(areaModulus*=3);
      
      v = bd.volume();
      a = bd.area(); 
      vred = 6.0*sqrt(M_PI)*v/std::pow(a,3.0/2.0);
      V = bd.prescribedVolume();
      A = bd.prescribedArea(); 
      Vred = 6.0*sqrt(M_PI)*V/std::pow(A,3.0/2.0);
      
      std::cout << "iter "<< iter << ":"<< std::endl
		<< "Volume = "<< bd.volume() << std::endl
		<< "Prescribed Volume = "<< bd.prescribedVolume() << std::endl
		<< "Area = "<< bd.area() << std::endl
		<< "Prescribed Area = "<< bd.prescribedArea() << std::endl
		<< "Reduced Volume = " << vred << std::endl
		<< "Prescribed Reduced Volume = " << Vred << std::endl
		<< "Energy = " << bd.energy() << std::endl
		<< "Stretching Energy = " << bd.stretchingEnergy() << std::endl
		<< "fixed pressure = " << bd.fixedPressure() << std::endl
		<< "fixed tension = " << bd.fixedTension() << std::endl;
      if( abs(vred-Vred) < 1.0e-3 && bd.stretchingEnergy() < 1.0e-6*bd.energy() )
	break;
    }
    sprintf(fname,"v-%f",Vred);
    model.print(fname);
  }
  return 0;
}

