#include <iostream>
#include <vector>
#include "SCElastic.h"
#include "StVenant.cc"
#include "ContractionWrapper.h"

int main(){

	std::cout << "\n \n \n";
	std::cout << "testing ContractionWrapper with StVenant Material ..." << std::endl;
	voom::StVenant stv(1.0, 1.0e9, 1.0/3.0);
	voom::ContractionWrapper<voom::StVenant> contraction(stv);
	contraction.setT(0.1);
	//
	// components of the deformation gradient tensor
	const double k11 = 1.6;
	const double k12 = 0.0;
	const double k13 = 0.0;
	const double k21 = 0.0;
	const double k22 = 1.0;
	const double k23 = 0.0;
	const double k31 = 0.0;
	const double k32 = 0.0;
	const double k33 = 1.0;
	tvmet::Matrix<double, 3, 3> F;
	F = k11, k12, k13, k21, k22, k23, k31, k32, k33;
	
	contraction.setDeformationGradient(F);
	std::cout << "deformation gradient tensor == " << contraction.deformationGradient() << std::endl;
	contraction.updateState(true, true, false);
	std::cout << "1st Piola stress = " << contraction.piolaStress() << std::endl;
	contraction.ConsistencyTest();
	
	return 0;
}
