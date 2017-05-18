// #include "C0Membrane.h"
// #include "C0MembraneShear.h"
#include "C0MembraneStretch.h"
#include "VoomMath.h"
#include "TriangleQuadrature.h"
#include "ShapeTri3.h"
// #include "EvansElastic_Skewed.h"
// #include "EvansElastic_SkewedMin.h"
#include "EvansElastic_Stretch.h"

using namespace voom;
using namespace std;

int main() {
  
  std::vector<Vector3D> X(3);
  X[0] = -4.04505, 0.754686, 20.0; // 0.0, 0.0, 0.0;
  X[1] = -3.03114, 1.004260, 20.0; // 2.0, 0.0, 0.0;
  X[2] =  3.13147, 0.119671, 20.0; // 0.0, 2.0, 0.0;

  std::vector<Vector3D> u(3);
  u[0] = 0.02146030, -0.00400383, -0.4233000;
  u[1] = 0.00967376, -0.00320540, -0.2549090;
  u[2] = 0.00962490, -0.000367821, -0.245511;

  srand(time(NULL));

  // typedef EvansElastic TestMaterial;
  // typedef EvansElastic_Skewed TestMaterial;
  // typedef C0Membrane<TriangleQuadrature, TestMaterial,  ShapeTri3> MemEl;

  // typedef EvansElastic_SkewedMin TestMaterial;
  // typedef C0MembraneShear MemEl;

  typedef EvansElastic_Stretch TestMaterial;  
  typedef C0MembraneStretch MemEl;
  MemEl::NodeContainer dnodes;

  int dof = 0, id = 0;
  for(id = 0; id < 3; id++)
  {
    Vector3D  x( X[id] );
    // perturb 
    x+=u[id];

    NodeBase::DofIndexMap idx(3);
    for(int j = 0; j < 3; j++) idx[j] = dof++;
    dnodes.push_back( new DeformationNode<3>(id,idx,X[id],x) );
  }
 
  TriangleQuadrature quad(1);
  Shape<2>::CoordinateArray s(0.0);
  ShapeTri3 shape(s);

  /*
  NodeBase::DofIndexMap idxM(1);
  idxM[0]=-1;
  MultiplierNode * pressureNode = new MultiplierNode(dnodes.size(), idxM, 0.0 );
  MultiplierNode * tensionNode = new MultiplierNode(dnodes.size(), idxM, 0.0 );
  */

  // Initialize stretch node (ScalarFieldNode type)
  id++; 
  NodeBase::DofIndexMap idE(1);   idE[0] = dof++;
  ScalarFieldNode<3>::PositionVector pE(0.0);
  ScalarFieldNode<3>* nodeEta = new ScalarFieldNode<3>(id, idE, pE, 1.006);
 
  // Initialize direction node (ScalarFieldNode type)
  id++;
  idE[0] = dof++;
  ScalarFieldNode<3>* nodeTheta = new ScalarFieldNode<3>(id, idE, pE, 0.1);

  // Initialize soft mode node (ScalarFieldNode type)
  id++;
  idE[0] = dof++;
  ScalarFieldNode<3>* nodePhi = new ScalarFieldNode<3>(id, idE, pE, 0.1);
 
  double mu = 1.2;
  double kS = 2.3;
  double angleOffset = 45.0;
  int WcType = 6;
  double WcConst[4] = {0.1, 0.005, 150.0, 1.06};
  TestMaterial stretch(mu, kS, nodeEta->point(), nodeTheta->point(), WcType, WcConst);
  // TestMaterial stretch(0.,0.,0.,mu,kS);
  // TestMaterial stretch(0.0, 0.0, 0.0, mu, kS, 1.0, 45.0, 0.0, 0.0);

  // MemEl* elem = new MemEl(quad, stretch, dnodes, pressureNode, tensionNode);
  MemEl * elem = new MemEl(dnodes, nodeEta, nodeTheta, nodePhi, angleOffset, quad, &stretch, &shape, true, true, true);
  // MemEl * elem = new MemEl(dnodes, nodeEta, quad, &stretch, &shape);


  elem->checkConsistency();
  cout << elem->strainEnergy() << endl;
  cout << elem->conformationalEnergy() << endl;

  delete elem;
  delete nodeEta;
  delete nodeTheta;
  for (int i = 0; i < dnodes.size(); i++) {delete dnodes[i];}

  return 0;

}

