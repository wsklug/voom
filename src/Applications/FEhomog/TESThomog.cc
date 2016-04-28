#include <string>
#include <iostream>
#include <vector>
#include "Body3D.h"
#include "Node.h"
#include "CompNeoHookean.h"
#include "StVenant.h"
#include "TetQuadrature.h"
#include "ShapeTet4CP.h"
#include <tvmet/Vector.h>
#include <fstream>


using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{

  // Create material
  double rho = 500.0;
  double E = 200.0;
  double nu = 0.4;
  // CompNeoHookean protein( rho, E, nu );
  StVenant protein1( rho, E, nu );
  StVenant protein2( rho, E, nu );
  vector<Material *> Protein;
  Protein.push_back(&protein1);
  Protein.push_back(&protein2);


  // Connectivity table
  vector<vector< int> > Connectivities(2, vector<int>(4, 0));
  Connectivities[0][0] = 0; Connectivities[0][1] = 1; Connectivities[0][2] = 2; Connectivities[0][3] = 3;
  Connectivities[1][0] = 0; Connectivities[1][1] = 1; Connectivities[1][2] = 2; Connectivities[1][3] = 4;



  // DefNodeContainer
  vector<DeformationNode<3> *> Nodes;
  Nodes.reserve(5);

  int dof = 0, j = 0, id = 0;
  NodeBase::DofIndexMap idx(3);
  DeformationNode<3>::Point x;

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 0.0; 
  DeformationNode<3>* a = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( a );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 1.0, 0.0; 
  DeformationNode<3>* b = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( b );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 1.0, 0.0, 0.0; 
  DeformationNode<3>* c = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( c );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0, 1.0; 
  DeformationNode<3>* d = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( d );

  id++;
  for(j = 0; j < 3; j++) idx[j]=dof++;
  x = 0.0, 0.0,-1.0; 
  DeformationNode<3>* e = new DeformationNode<3>(id, idx, x);
  Nodes.push_back( e );

  
  // Assign linear deformation gradient
  Tensor3D F(0.0), C(0.0), Csquare(0.0);
  F(0,0) = 1.0; F(0,1) = 2.0; F(0,2) =-1.0;
  F(1,0) =-2.0; F(1,1) = 1.5; F(1,2) = 5.0;
  F(2,0) = 1.2; F(2,1) = 3.0; F(2,2) = 1.0;
 
  for (int i=0; i< Nodes.size(); i++)
  {
    tvmet::Vector<double,3> PointTaken, PointSet(0.0);

    PointTaken = Nodes[i]->point();

    PointSet(0) = F(0,0)*PointTaken(0) + F(0,1)*PointTaken(1) + F(0,2)*PointTaken(2);
    PointSet(1) = F(1,0)*PointTaken(0) + F(1,1)*PointTaken(1) + F(1,2)*PointTaken(2);
    PointSet(2) = F(2,0)*PointTaken(0) + F(2,1)*PointTaken(1) + F(2,2)*PointTaken(2);
      
    Nodes[i]->setPoint(PointSet);
  }
    
      // Compute invariants
      double detF = determinant(F); 
      double detF_TwoThird = pow(detF, -2.0/3.0);
      int i = 0, k = 0;
      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  for (k = 0; k < 3; k++) {
	    C(i,j) += F(k,i)*F(k,j);
	  }
	}
      }
      for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	  for (k = 0; k < 3; k++) {
	    Csquare(i,j) += C(i,k)*C(k,j);
	  }
	}
      }

      double trC = C(0,0) + C(1,1) + C(2,2);

      cout << "Expected inv1 = " << trC *detF_TwoThird << endl;
      cout << "Expected inv2 = " << 0.5*(pow(trC,2.0) - Csquare(0,0) - Csquare(1,1) - Csquare(2,2))*pow(detF_TwoThird, 2.0) << endl;
      cout << "Expected inv3 = " << detF << endl; //*detF << endl;

      


  // Create quadrature rule and shape function
  TetQuadrature Quad(1);
  ShapeTet4 Sh;
 
  
  
  // Create body
  Body3D BodyTest(Protein, Connectivities, Nodes, Quad, Sh, 1.0);
  cout << "I initialized the body" << endl;
  
  // Should do nothing - nothing to be reset - reference nodal positions did not change
  BodyTest.reset();

  // d->setPoint(2,1.1);
  BodyTest.checkConsistency();

  // Should do nothing - nothing to be reset - reference nodal positions did not change
  BodyTest.reset();

  int f = 0;
  vector<pair<Vector3D, vector<double > > > Invariants = BodyTest.invariants(f);
  for (int i = 0; i<Invariants.size(); i++) {
    cout << (Invariants[i]).first << "  ";
    for (int j =0; j<3; j++) {
      cout << ((Invariants[i]).second)[j] << " "; 
    }
    cout << endl;
  }
  cout << endl;
  cout << "Element with detF < 0.0 = " << f << endl;
  
  for (int k = 0; k<Nodes.size(); k++)
  {
    delete Nodes[k];
  }
	

  return 0;
}



