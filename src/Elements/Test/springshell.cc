#include "SpringShell.h"

using namespace voom;
using std::cout;
using std::setw;
using std::endl;

int main() {

  //
  // 2-------1
  // | (0) / |
  // |   /   |
  // | / (1) |
  // 0-------3
  //

  std::vector<Vector3D> X(3);
  X[0] = 0.0, 0.0, 0.0;
  X[1] = 2.0, 2.0, 0.0;
  X[2] = 0.0, 2.0, 1.0;
  X[3] = 2.0, 0.0, 1.0;

  const int nNodes = 4;

  srand(time(0));

  SpringShell::DefNodeContainer dnodes;

  int dof=0;
  for(int id=0; id<nNodes; id++) {
    //X[id] *= 0.01;

    Vector3D  x( X[id] );

    // perturb current positions by random displacement
    for(int i=0; i<3; i++) x(i) += 0.01*(rand()/((double)RAND_MAX)-0.5);

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

  double scalingFactor=1.0e0;

  double k = 1.0;//4.0e4;
  double theta0 = 0.1;

  SpringShell * elem = new SpringShell(dnodes, k, theta0);

  bool passed = elem->checkConsistency();

  cout << "Checking consistency: " << passed << endl;

  delete elem;
  
  for(int id=0; id<dnodes.size(); id++) {
    delete dnodes[id];
  }

  cout << "Bye bye." << endl;
  return 0;

}
