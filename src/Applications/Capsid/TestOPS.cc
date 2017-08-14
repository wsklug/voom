/*
 * OPSAsphericity.cc
 *
 *  Created on: Aug 7, 2017
 *      Author: amit
 */
#include <string>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <vector>
#include <fstream>

#include <tvmet/Vector.h>
#include <limits>
#include "Node.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "OPSBody.h"

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>
#include "HelperFunctions.h"

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[])
{
  clock_t t1, t2, t3;
  t1 = clock();
  if (argc != 2) {
    //cout << "usage: " << argv[0] << " <filename>\n";
    //return -1;
  }

  int dof = 0;
  Vector3D x1(0.0), v1(0.0), v2(0.0), x2(0.0);
  NodeBase::DofIndexMap idx1(6), idx2(6);
  OPSNode *n1, *n2;
  vector<OPSNode*> nodes;

  x1 = 0.5, 0.75, 0.85;
  v1 = 0.3, 0.4, 0.8660254037844386;

  x2 = 1.5, 1.75, 1.85;
  v2 = 0.4, 0.8660254037844386, 0.3;

  for (int j = 0; j < 6; j++) idx1[j] = dof++;
  for (int j = 0; j < 6; j++) idx2[j] = dof++;

  n1 = new OPSNode(0, idx1, x1);
  n2 = new OPSNode(1, idx2, x2);
  n1->setDeformedRotationVector(v1);
  n2->setDeformedRotationVector(v2);

  nodes.push_back(n1);
  nodes.push_back(n2);

  struct OPSParams prop = {};
  OPSBody bd(nodes, prop, 5.0);

  double morseEn, psiVal, phi_pVal, phi_nVal, phi_cVal;
  Vector3D dMdXi, dMdXj, dSdXi, dSdXj, dSdVi, dPdVi, dPdXi, dPdXj;
  Vector3D dNdVi, dNdVj, dCdVi, dCdVj, dCdXi, dCdXj;

  morseEn = bd.morse(x1, x2);
  dMdXi = bd.DmorseDxi(x1,x2);
  dMdXj = bd.DmorseDxj(x1,x2);

  psiVal = bd.psi(v1,x1,x2);
  dSdVi = bd.DpsiDvi(v1,x1,x2);
  dSdXi = bd.DpsiDxi(v1,x1,x2);
  dSdXj = bd.DpsiDxj(v1,x1,x2);

  phi_pVal = bd.phi_p(v1,x1,x2);
  dPdVi = bd.Dphi_pDvi(v1,x1,x2);
  dPdXi = bd.Dphi_pDxi(v1,x1,x2);
  dPdXj = bd.Dphi_pDxj(v1,x1,x2);

  phi_nVal = bd.phi_n(v1,v2);
  dNdVi = bd.Dphi_nDvi(v1,v2);
  dNdVj = bd.Dphi_nDvj(v1,v2);

  phi_cVal = bd.phi_c(v1,v2,x1,x2);
  dCdVi = bd.Dphi_cDvi(v1,v2,x1,x2);
  dCdVj = bd.Dphi_cDvj(v1,v2,x1,x2);
  dCdXi = bd.Dphi_cDxi(v1,v2,x1,x2);
  dCdXj = bd.Dphi_cDxj(v1,v2,x1,x2);

  cout<<"Printing the values of potential functions and derivatives: " << endl;
  cout<<"\t Morse = " << morseEn << endl;
  cout<<"\t Psi = " << psiVal << endl;
  cout<<"\t Phi_p = "<< phi_pVal << endl;
  cout<<"\t Phi_n = "<< phi_nVal << endl;
  cout<<"\t Phi_c = "<< phi_cVal << endl;
  cout<<"\t dMdXi = " << dMdXi[0] << "," << dMdXi[1] << "," << dMdXi[2] << endl;
  cout<<"\t dMdXj = " << dMdXj[0] << "," << dMdXj[1] << "," << dMdXj[2] << endl;
  cout<<"\t dSdVi = " << dSdVi[0] << "," << dSdVi[1] << "," << dSdVi[2] << endl;
  cout<<"\t dSdXi = " << dSdXi[0] << "," << dSdXi[1] << "," << dSdXi[2] << endl;
  cout<<"\t dSdXj = " << dSdXj[0] << "," << dSdXj[1] << "," << dSdXj[2] << endl;
  cout<<"\t dPdVi = " << dPdVi[0] << "," << dPdVi[1] << "," << dPdVi[2] << endl;
  cout<<"\t dPdXi = " << dPdXi[0] << "," << dPdXi[1] << "," << dPdXi[2] << endl;
  cout<<"\t dPdXj = " << dPdXj[0] << "," << dPdXj[1] << "," << dPdXj[2] << endl;
  cout<<"\t dNdVi = " << dNdVi[0] << "," << dNdVi[1] << "," << dNdVi[2] << endl;
  cout<<"\t dNdVj = " << dNdVj[0] << "," << dNdVj[1] << "," << dNdVj[2] << endl;
  cout<<"\t dCdVi = " << dCdVi[0] << "," << dCdVi[1] << "," << dCdVi[2] << endl;
  cout<<"\t dCdVj = " << dCdVj[0] << "," << dCdVj[1] << "," << dCdVj[2] << endl;
  cout<<"\t dCdXi = " << dCdXi[0] << "," << dCdXi[1] << "," << dCdXi[2] << endl;
  cout<<"\t dCdXj = " << dCdXj[0] << "," << dCdXj[1] << "," << dCdXj[2] << endl;

  delete n1, n2;
  t3 = clock();
  float diff = ((float)t3 - (float)t2);
  std::cout << "Post-processing execution time: " << diff / CLOCKS_PER_SEC
	    << " seconds" << std::endl;

}
