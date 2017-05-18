#include "VoomMath.h"
#include "PotentialElement.h"
#include "LennardJones.h"
#include "SpringPotential.h"

using namespace voom;
using namespace std;

int main() {
  
  // Declare nodes
  std::vector<Vector3D> X(3);
  X[0] = -4.04505, 0.754686, 20.0; // 0.0, 0.0, 0.0;
  X[1] = -3.03114, 1.004260, 20.0; // 2.0, 0.0, 0.0;
  X[2] =  3.13147, 0.119671, 20.0; // 0.0, 2.0, 0.0;
  X[3] =  0.0, 0.0, 0.0;

  std::vector<Vector3D> u(3);
  u[0] = 0.02146030, -0.00400383, -0.4233000;
  u[1] = 0.00967376, -0.00320540, -0.2549090;
  u[2] = 0.00962490, -0.000367821, -0.245511;

  set<DeformationNode<3> * > domain;
  int dof = 0, id = 0;
  for(id = 0; id < 3; id++)
  {
    Vector3D  x( X[id] );
    x+=u[id];
   
    NodeBase::DofIndexMap idx(3);
    for(int j = 0; j < 3; j++) idx[j] = dof++;
    domain.insert( new DeformationNode<3>(id,idx,X[id],x) );
  }

  NodeBase::DofIndexMap idx(3);
  for(int j = 0; j < 3; j++) idx[j] = dof++;
  DeformationNode<3> * center = new DeformationNode<3>(3,idx,X[3]);



  // Material
  srand(time(NULL));
  double epsilon = double(rand())/double(RAND_MAX);
  double sigma = double(rand())/double(RAND_MAX);
  
  LennardJones Mat1(sigma, epsilon);
  SpringPotential Mat2(sigma, epsilon);


  // Element
  PotentialElement MyEl1(&Mat1, center, domain);
  MyEl1.compute(true, false, false);
  cout << "EnergyEl1 = " << MyEl1.energy() << endl;
  MyEl1.checkConsistency();

  PotentialElement MyEl2(&Mat2, center, domain);
  MyEl2.compute(true, false, false);
  cout << "EnergyEl2 = " << MyEl2.energy() << endl;
  MyEl2.checkConsistency();

  delete center;
  for (set<DeformationNode<3> *>::iterator pNode = domain.begin();
        pNode != domain.end(); pNode++)
  {
    delete *pNode;
  }
}
