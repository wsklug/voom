#include <vector>
#include <iostream>
#include "SCElastic.h"

int main()
{
  VOOM::SCElastic bilayer( 1.0, 0.0 );
  tvmet::Vector<Vector3D,2> a;
  tvmet::Matrix<Vector3D,2,2> aPartials;

  a(0) = 1.0, 0.0, 0.0; 
  a(1) = 0.0, 1.0, 0.0; 

  aPartials(0,0) = 0.0, 0.0, -1.0;
  aPartials(1,1) = 0.0, 0.0, 1.0;
  aPartials(0,1) = 0.0, 0.0, 0.0;
  aPartials(1,0) = 0.0, 0.0, 0.0;

  bilayer.setGeometry( VOOM::ShellGeometry( a, aPartials ) );
  bilayer.updateState(true,true,false);

  for(int i=0; i<3; i++)
    std::cout << "n[" << i << "]=" << bilayer.stressResultants()(i) 
	      << std::endl;
  for(int i=0; i<2; i++)
    std::cout << "m[" << i << "]=" << bilayer.momentResultants()(i) 
	      << std::endl;

  std::cout << "H=" << bilayer.meanCurvature() << std::endl;
  std::cout << "K=" << bilayer.gaussianCurvature() << std::endl;
  std::cout << "E=" << bilayer.energy() << std::endl;
  return 0;
}
