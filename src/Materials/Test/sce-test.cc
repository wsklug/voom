// -*- C++ -*-
#include <vector>
#include <iostream>
#include "SCElastic.h"

int main() {

  std::cout << "testing SCElastic material..." << std::endl;
	
  voom::SCElastic bilayer( 1.0, 1.0 );
  tvmet::Vector<Vector3D,2> a;
  tvmet::Matrix<Vector3D,2,2> aPartials;

  a(0) = 1.0, 0.0, 0.10; 
  a(1) = 0.0, 1.0, 0.10; 

  aPartials(0,0) = 0.30, -0.40, -1.0;
  aPartials(1,1) = 0.10, 0.20, 1.50;
  aPartials(0,1) = 0.20, 0.10, 0.40;
  aPartials(1,0) = aPartials(0,1);
  voom::ShellGeometry geometry(a, aPartials);
  bilayer.setGeometry( geometry );
  bilayer.updateState(true,true,false);

  for(int i=0; i<3; i++)
    std::cout << "n[" << i << "]=" << bilayer.stressResultants()(i) 
	      << std::endl;
  for(int i=0; i<2; i++)
    std::cout << "m[" << i << "]=" << bilayer.momentResultants()(i) 
	      << std::endl;

  std::cout << "H=" << bilayer.meanCurvature() << std::endl;
  std::cout << "K=" << bilayer.gaussianCurvature() << std::endl;
  std::cout << "W=" << bilayer.energyDensity() << std::endl;

  double h=1.0e-8;
  tvmet::Vector<Vector3D,2> n(Vector3D(0.0));
  tvmet::Vector<Vector3D,2> m(Vector3D(0.0));
  for(int alpha=0; alpha<2; alpha++) {
    for(int i=0; i<3; i++) {
      a(alpha)(i) += h;
      bilayer.setGeometry( voom::ShellGeometry(a,aPartials) );
      bilayer.updateState(true,false,false);
      n(alpha)(i) = bilayer.energyDensity();
      a(alpha)(i) -= 2.0*h;
      bilayer.setGeometry( voom::ShellGeometry(a,aPartials) );
      bilayer.updateState(true,false,false);
      n(alpha)(i) -= bilayer.energyDensity();
      a(alpha)(i) += h;

      n(alpha)(i) /= 2.0*h;

      bilayer.setGeometry ( geometry );
      bilayer.updateState(true,false,false);

      geometry.dPartials()(alpha)(i) += h;
      bilayer.setGeometry( geometry );
      bilayer.updateState(true,false,false);
      m(alpha)(i) = bilayer.energyDensity();
      geometry.dPartials()(alpha)(i) -= 2.0*h;
      bilayer.setGeometry( geometry );
      bilayer.updateState(true,false,false);
      m(alpha)(i) -= bilayer.energyDensity();
      geometry.a()(alpha)(i) += h;
      
      m(alpha)(i) /= 2.0*h;

      bilayer.setGeometry ( geometry );
      bilayer.updateState(true,false,false);

    }
  }

  bilayer.setGeometry( geometry );
  bilayer.updateState(false,true,false);
  double normStress=0;
  double normMoment=0;
  double errorStress=0;
  double errorMoment=0;
  for(int alpha=0; alpha<2; alpha++) {
    for(int i=0; i<3; i++) {
      normStress = std::max( normStress, std::abs(n(alpha)(i)));
      normMoment = std::max( normMoment, std::abs(m(alpha)(i)));
      double tmp = n(alpha)(i) - bilayer.stressResultants()(alpha)(i);
      errorStress += tmp*tmp;
      tmp = m(alpha)(i) - bilayer.momentResultants()(alpha)(i);
      errorMoment += tmp*tmp;
    }
  }
  errorStress = sqrt(errorStress);
  errorMoment = sqrt(errorMoment);

  std::cout << "errorStress = " << errorStress << '\t'
	    << "h*normStress = " << h*normStress << '\t' << std::endl
	    << "errorMoment = " << errorMoment << '\t'
	    << "h*normMoment = " << h*normMoment << '\t' << std::endl;

  for(int i=0; i<2; i++)
    std::cout << "n[" << i << "]=" << n(i) 
	      << std::endl;
  for(int i=0; i<2; i++)
    std::cout << "m[" << i << "]=" << m(i) 
	      << std::endl;
  return 0;
}
