#include "C0Membrane.h"
#include "VoomMath.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
#include "FVK.h"
#include "Graphene.h"

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
//   u[0] = 0.0214603, -0.00400383, -0.4233;
//   u[1] = 0.00967376, -0.0032054, -0.254909;
//   u[2] = 0.0096249, -0.000367821, -0.245511;
  u[0] = 0.0,0.0,0.0;
  u[1] = 0.0,0.0,0.0;
  u[2] = 0.0,0.0,0.0;


  srand(time(0));

  // typedef C0Membrane<TriangleQuadrature,FVK,ShapeTri3> ELTYPE;
  typedef C0Membrane<TriangleQuadrature,Graphene,ShapeTri3> ELTYPE;
  ELTYPE::NodeContainer dnodes;

  int dof=0;
  for(int id=0; id<3; id++) {
    //X[id] *= 0.01;

    Vector3D  x( X[id] );

    // perturb current positions by random displacement
    for(int i=0; i<3; i++) x(i) += 0.01*(rand()/((double)RAND_MAX)-0.5);
    // x+=u[id];

    cout << id << " X     x" <<  endl 
	 << X[id](0) << ',' << x(0) << endl
	 << X[id](1) << ',' << x(1) << endl
	 << X[id](2) << ',' << x(2) << endl
	 << endl;

    NodeBase::DofIndexMap idx(3);
    for(int j=0; j<3; j++) idx[j]=dof++;
    dnodes.push_back( new DeformationNode<3>(id,idx,X[id],x) );

  }

  double scalingFactor=1.0e0;

  double Y = 100.0;//4.0e4;
  double KC = 0.0;
  double nu = 0.3;//1.0/3.0;
  
  double KG =-KC ;
  double C0 = 0.0;

  double c111=-2724.7;
  double c222=-2523.2;
  double c112=-591.1;

  Y *= scalingFactor;

  // FVK mat( KC, KG, C0, Y, nu );	
  Graphene mat( Y, nu, c111, c222, c112 );	
  
  unsigned int quadOrder = 1;

  TriangleQuadrature quad(1);

  ShapeTri3 shape( Shape<2>::CoordinateArray(0.0) );

  ELTYPE* elem = new ELTYPE(quad, mat, dnodes, 0, 0);

  elem->checkConsistency();

  cout << "Bye bye." << endl;
  return 0;

}
