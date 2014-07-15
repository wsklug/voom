#include "C0MembraneGL.h"
#include "VoomMath.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "MembraneGLImplicitMass.h"

using namespace voom;
using std::cout;
using std::setw;
using std::endl;

int main() {
  
  std::vector<Vector3D> X(3);
  X[0] = -4.04505, 0.754686, 20.0;//0.0, 0.0, 0.0;
  X[1] = -3.03114, 1.00426, 20.0;//2.0, 0.0, 0.0;
  X[2] = 3.13147, 0.119671, 20.0;//0.0, 2.0, 0.0;

  std::vector<Vector3D> u(3);
  u[0] = 0.0214603, -0.00400383, -0.4233;
  u[1] = 0.00967376, -0.0032054, -0.254909;
  u[2] = 0.0096249, -0.000367821, -0.245511;

  srand(time(0));

  C0MembraneGL::DefNodeContainer dnodes;
  C0MembraneGL::GLNodeContainer gnodes;

  int dof=0;
  for(int id=0; id<3; id++) {
    //X[id] *= 0.01;

    Vector3D  x( X[id] );

    // perturb current positions by random displacement
    for(int i=0; i<3; i++) x(i) += 0.01*(rand()/((double)RAND_MAX)-0.5);
    x+=u[id];

    cout << id << " (" 
	 << X[id](0) << ','
	 << X[id](1) << ','
	 << X[id](2) << ") ("
	 << x(0) << ','
	 << x(1) << ','
	 << x(2) << ")"
	 << endl;

    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    dnodes.push_back( new DeformationNode<3>(id,idx,X[id],x) );

  }

  for(int id=0; id<3; id++) {  
    double rho = 1.0;//rand()/((double)RAND_MAX);
    cout << id << " " << rho << endl;
    NodeBase::DofIndexMap idx(1);
    idx[0]=dof++;
    gnodes.push_back( new ScalarFieldNode<3>(id,idx,X[id],rho) );
  }

  double scalingFactor=1.0e0;

  double Y = 100.0;//4.0e4;
  double KC = 0.0;
  double nu = 0.3;//1.0/3.0;
  
  double KG =-KC ;
  double C0 = 0.0;

  double dg = 0.0;
  double g0= 1.5;
  double Gamma=1.5;
  double muA=1.0;
  int formulation=0;

  Y *= scalingFactor;
  dg *= scalingFactor;
  g0 *= scalingFactor;
  Gamma *= scalingFactor;
  muA *= scalingFactor;

  GLElastic mat( KC, KG, C0, Y, nu, Gamma, g0, dg, muA, formulation );
  
  unsigned int quadOrder = 1;

  TriangleQuadrature quad(1);

  ShapeTri3 shape( Shape<2>::CoordinateArray(0.0) );

  C0MembraneGL* elem = new C0MembraneGL(dnodes, gnodes, &mat, &quad, &shape);

  elem->checkConsistency();

  double Mx0=1.0e0, Mx1=0.0, Mrho0=1.0e0;
  
  Mx0 *= scalingFactor;
  Mx1 *= scalingFactor;
  Mrho0 *= scalingFactor;

  TriangleQuadrature quad2(2);
  
  MembraneGLImplicitMass * mass 
    = new MembraneGLImplicitMass(dnodes, gnodes, &quad2, &shape, 
				 Mx0, Mx1, Mrho0);
  double dt=1.0e-4;
  mass->step(dt);

  // perturb current positions by random displacement
  for(int a=0; a<3; a++) {
    for(int i=0; i<3; i++) {
      dnodes[a]->addPoint(i, 0.01*(rand()/((double)RAND_MAX)-0.5) );
    }
    gnodes[a]->setPoint(rand()/((double)RAND_MAX));
  }
  mass->checkConsistency();
  cout << "Bye bye." << endl;
  return 0;

}
