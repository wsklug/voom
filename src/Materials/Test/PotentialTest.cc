#include "LennardJones.h"
#include "LennardJonesFT.h"
#include "SpringPotential.h"
#include "SpringPotentialSQ.h"
#include "Morse.h"
#include "ProteinMorse.h"
#include "NodeBase.h"

using namespace voom;

int main()  
{
  srand (time(NULL));
  double epsilon = double(rand())/double(RAND_MAX);
  double sigma = double(rand())/double(RAND_MAX);
  double Rshift = double(rand())/double(RAND_MAX);

  int id = 0, dof = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;

  for(uint j = 0; j < 3; j++) idx[j] = dof++;
  x(0) = 0.3; x(1) = 0.1; x(2) = 0.2;
  DeformationNode<3>* nodeA = new DeformationNode<3>(id,idx,x);
  cout << "ForceA = " << nodeA->force() << endl;

  id = 1;
  for(uint j = 0; j < 3; j++) idx[j] = dof++;
  x(0) = 1.3; x(1) = 1.1; x(2) = 1.2;
  DeformationNode<3>* nodeB = new DeformationNode<3>(id,idx,x);
  cout << "ForceB = " << nodeB->force() << endl;

  {
    LennardJones LennPotential(epsilon, sigma);
    LennPotential.ConsistencyTest(nodeA, nodeB, 1.0e-7, 1.0e-7);
    LennPotential.updateState(nodeA, nodeB, true,true,false);
    cout << "LennerdJones energy = " << LennPotential.energy() << endl;
    cout << endl << " ---------------------------------------------- " << endl;
  }
 
  {
    LennardJonesFT LennPotential42(epsilon, sigma, Rshift);
    LennPotential42.ConsistencyTest(nodeA, nodeB, 1.0e-7, 1.0e-7);
    LennPotential42.updateState(nodeA, nodeB, true,true,false);
    cout << "LennerdJones42 energy = " << LennPotential42.energy() << endl;
    cout << endl << " ---------------------------------------------- " << endl;
  }

  {
    SpringPotential SpPotential(epsilon, sigma);
    SpPotential.ConsistencyTest(nodeA, nodeB, 1.0e-7, 1.0e-7);
    SpPotential.updateState(nodeA, nodeB, true,true,false);
    cout << "Spring energy = " << SpPotential.energy() << endl;
    cout << endl << " ---------------------------------------------- " << endl;
  }

  {
    SpringPotentialSQ SpPotentialSQ(epsilon, sigma);
    SpPotentialSQ.ConsistencyTest(nodeA, nodeB, 1.0e-7, 1.0e-7);
    SpPotentialSQ.updateState(nodeA, nodeB, true,true,false);
    cout << "Spring energy = " << SpPotentialSQ.energy() << endl;
    cout << endl << " ---------------------------------------------- " << endl;
  }

  {
    Morse MorsePotential(epsilon, sigma, Rshift); 
    MorsePotential.ConsistencyTest(nodeA, nodeB, 1.0e-7, 1.0e-7);
    MorsePotential.updateState(nodeA, nodeB, true,true,false);
    cout << "Morse energy = " << MorsePotential.energy() << endl;
    cout << endl << " ---------------------------------------------- " << endl;
  }

  {
    ProteinMorse PMorse(epsilon, sigma, Rshift); 
    ProteinNode pA(nodeA), pB(nodeB);
    cout << "Morse energy = " << PMorse.computeEnergy(&pA, &pB) << endl;
    cout << endl << " --------------------------------------- " << endl;
  }

  delete nodeA;
  delete nodeB;
  
  return 0;
}
