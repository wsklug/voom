#include "../LMEshape.h"
#include "../LMEshape2.h"
#include <iostream>
//#include <rand>

using namespace std;
using namespace voom;

int main()
{
  srand(time(0));

  double beta = 0.8;
  double searchR = 2.0;
  double tol = 1.0e-12;
  double nItMax = 20;

  cout << "Test LME in 3D" << endl;

  {
    LMEshape LME3D(beta, searchR, tol, nItMax);

    LMEshape::CoordinateArray s;
    s(0) = 0.05;
    s(1) = 0.55;
    s(2) = 0.23;



    vector<DeformationNode<3>* > Nodes;
    vector<int> Dof(3,1);
    tvmet::Vector<double,3> p(0.0);
    
    tvmet::Vector<double,3> XA(0.0, 0.0, 0.0);
    DeformationNode<3> A(0, Dof, XA, p);
    Nodes.push_back(&A);

    tvmet::Vector<double,3> XB(1.0, 0.0, 0.0);
    DeformationNode<3> B(1, Dof, XB, p);
    Nodes.push_back(&B);

    tvmet::Vector<double,3> XC(1.0, 1.0, 0.0);
    DeformationNode<3> C(2, Dof, XC, p);
    Nodes.push_back(&C);

    tvmet::Vector<double,3> XD(0.0, 1.0, 0.0);
    DeformationNode<3> D(3, Dof, XD, p);
    Nodes.push_back(&D);

    tvmet::Vector<double,3> XE(0.0, 0.0, 1.0);
    DeformationNode<3> E(4, Dof, XE, p);
    Nodes.push_back(&E);

    tvmet::Vector<double,3> XF(1.0, 0.0, 1.0);
    DeformationNode<3> F(5, Dof, XF, p);
    Nodes.push_back(&F);

    tvmet::Vector<double,3> XG(1.0, 1.0, 1.0);
    DeformationNode<3> G(6, Dof, XG, p);
    Nodes.push_back(&G);

    tvmet::Vector<double,3> XH(0.0, 1.0, 1.0);
    DeformationNode<3> H(7, Dof, XH, p);
    Nodes.push_back(&H);
  
    LME3D.compute(s, Nodes, true, true);
    
    double Sum = LME3D.checkPartitionOfUnity();
  
    cout << "The sum of the shape functions over all nodes is = " << Sum << endl;

    vector<double > Na  = LME3D.functions();
    vector<double > DNx  = LME3D.xderivative();
    vector<double > DNy  = LME3D.yderivative();
    vector<double > DNz  = LME3D.zderivative();

    for (int i = 0; i < Na.size(); i++)
      { 
	cout << "i = " << i << endl;
	cout << Na[i] << endl;
	cout << DNx[i] << "  " << DNy[i] << "  " << DNz[i] << endl;
      }
  }



  cout << "Test LME in 2D" << endl;

  {
    LMEshape2 LME2D(beta, searchR, tol, nItMax);

    LMEshape2::CoordinateArray s;
    s(0) = 0.0;
    s(1) = 0.0;



    vector<DeformationNode<2>* > Nodes;
    vector<int> Dof(2,1);
    tvmet::Vector<double,2> p(0.0);
    
    tvmet::Vector<double,2> XA(0.0, 0.0);
    DeformationNode<2> A(0, Dof, XA, p);
    Nodes.push_back(&A);

    tvmet::Vector<double,2> XB(1.0, 0.0);
    DeformationNode<2> B(1, Dof, XB, p);
    Nodes.push_back(&B);

    tvmet::Vector<double,2> XC(1.0, 1.0);
    DeformationNode<2> C(2, Dof, XC, p);
    Nodes.push_back(&C);

    tvmet::Vector<double,2> XD(0.0, 1.0);
    DeformationNode<2> D(3, Dof, XD, p);
    Nodes.push_back(&D);

    tvmet::Vector<double,2> XE(0.5, 0.5);
    DeformationNode<2> E(4, Dof, XE, p);
    Nodes.push_back(&E);
  
    LME2D.compute(s, Nodes, true);
    
    double Sum = LME2D.checkPartitionOfUnity();
  
    cout << "The sum of the shape functions over all nodes is = " << Sum << endl;

    LMEshape2::FunctionContainer N =  LME2D.functions();
    LMEshape2::DerivativeContainer DN =  LME2D.derivatives();

    for (int i = 0; i < N.size(); i++)
      { 
	cout << "i = " << i << endl;
	cout << N[i] << endl;
	cout << DN[i][0] << "  " << DN[i][1] << endl;
      }
  }




  
  return 0;
}
