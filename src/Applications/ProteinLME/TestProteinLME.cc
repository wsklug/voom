// This file reads in the configurations from MD
// either in .dat or .pdf formats
// Then calculates the deformation gradient at each node
// Then calculates the three invariants at each node
// and writes them to the file invariants1.dat *2.dat and *3.dat
// in the format : at each configuration a row of invariants at all the nodes
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace std;
#include "voom.h"
#include "LMEbodyQP.h"
#include "LMEshape.h"
#include "NodeBase.h"
#include "Node.h"
#include "Model.h"
#include "CompNeoHookean.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>

using namespace voom;

int main(int argc, char* argv[])
{ 
  bool verbose=true;
  
  // Input parameters
  ifstream ifs;
  double searchR = 2.0;
  double beta = 0.8;
  int nItMax = 100;
  double tol = 1.0e-10;

    // Create vector of nodes
    int dof = 0;
    std::vector< DeformationNode<3>* > Nodes;
    std::vector< DeformationNode<3>* > QP;
    
    double Lx = 1.5, Ly = 2.0, Lz = 1.0;
    int nx = 6, ny = 6, nz = 5;
    double DeltaX = Lx/double(nx-1), DeltaY = Ly/double(ny-1), DeltaZ = Lz/double(nz-1);
    int npts = nx*ny*nz;
    Nodes.reserve(npts);
    vector<double > supp_size(npts, 1.0);
    vector<double > volume(npts, 1.0);

    NodeBase::DofIndexMap idx(3);
    DeformationNode<3>::Point x;
    DeformationNode<3>::Point xqp;
    xqp(0) = 0.0;
    xqp(1) = 0.0;
    xqp(2) = 0.0;

    // Generate points
    int id = 0;
    for(int i = 0; i < nx; i++) {
      for(int j = 0; j < ny; j++) {
	for(int k = 0; k < nz; k++) {
	  id = i;
	  
	  x(0) = i*DeltaX+0.001;//*i-0.004*j*j;
	  x(1) = j*DeltaY+0.001*k*k*k;
	  x(2) = k*DeltaZ-0.002*i*j*k;
	  
	  for(int m=0; m<3; m++) idx[m] = dof++;
	  DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);

          xqp(0) += x(0);
	  xqp(1) += x(1);
	  xqp(2) += x(2);

	  Nodes.push_back( n );
	}
      }
    }
    double Gsupp = (DeltaX + DeltaY + DeltaZ)/3.0;



    for(int j=0; j<3; j++) idx[j] = dof++;
    xqp(0) /= npts;   xqp(1) /= npts;   xqp(2) /= npts;

    DeformationNode<3>* Xqp = new DeformationNode<3>(id++, idx, xqp);
    QP.reserve(1);
    QP.push_back(Xqp);
    cout << "QP->position() = " << QP[0]->position() << endl;

    beta /= pow(Gsupp,2.0);
    cout << "Number of nodes: " << Nodes.size() << endl;
    cout << "Avg supp size: " << Gsupp << endl;
    cout << "LME const beta used: " << beta << endl;
 



  // Create material, values don't matter here
  double rho = 1.0;
  double E = 200.0;
  double nu = 0.4;
   
  // Create and initialize body
  typedef CompNeoHookean MaterialType;
  MaterialType protein( rho, E, nu );
  
  // create Body
  typedef LMEbodyQP<MaterialType, LMEshape> body;
  body LMEbd(protein, Nodes, QP, volume, supp_size, beta, searchR, tol, nItMax);
 
  // LMEbd.setOutput(paraview);
  // LMEbd.compute(true,true,false);
  
  std::cout << std::endl;
  std::cout << "Body energy = " << LMEbd.energy() << std::endl;
  std::cout << "Body volume = " << LMEbd.volume() << std::endl;
  std::cout << std::endl;






  // Assign linear deformation gradient
  Tensor3D F(0.0), C(0.0), Csquare(0.0);
  F(0,0) = 1.0; F(0,1) = 2.0; F(0,2) =-1.0;
  F(1,0) =-2.0; F(1,1) = 1.5; F(1,2) = 5.0;
  F(2,0) = 1.2; F(2,1) = 3.0; F(2,2) = 1.0;
 
  for (int i=0; i< npts; i++)
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
      int i = 0, j = 0, k = 0;
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





    // Compute the invariants
    std::vector<double> I1(QP.size(),0.), I2(QP.size(),0.), I3(QP.size(),0.);
    // std::cout<< "Calculating the invariants" << std::endl;
    LMEbd.cal_invariants(I1,I2,I3);
 
    
    
    for(int iqp = 0; iqp < QP.size(); iqp++)
    {
      cout << QP[iqp]->position() << endl;
      cout << I1[iqp] << " " <<  I2[iqp] << " " <<  I3[iqp] << " " << endl;
    }
    
  
  std::cout<<"Everything done!"<<std::endl;
  return 0;

}

