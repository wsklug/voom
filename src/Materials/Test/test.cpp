#include "Hookean2D.h"
#include "Hookean.h"
#include "LinearizedMaterial-impl.h"
#include "DefinedTypes.h"
#include <iostream>

int main()
{
  voom::Hookean2D planeStrain(voom::Hookean2D::planeStrain, 10.0, 0.3);
  planeStrain.checkConsistency();
  voom::Hookean2D planeStress(voom::Hookean2D::planeStress, 10.0, 0.3);
  planeStress.checkConsistency();
  voom::Hookean threeDimensional(10.0, 0.3);
  threeDimensional.checkConsistency();

  planeStrain.checkIsotropy();
  planeStress.checkIsotropy();
  
  voom::Hookean mat(27.0, 0.3);
  voom::Hookean::VoigtVector e;
  e = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  mat.setStrain(e);
  mat.updateState(true,true,false);
  std::cout<<"sigma = "<<mat.stress()<<std::endl
	   <<"w = "<<mat.energyDensity()<<std::endl;

  return 0;
}


// #include <vector>
// #include <iostream>
// #include "SCElastic.h"
// #include "StVenant.h"


// 	std::cout << "\n \n\n" ;
// 	std::cout << "testing StVenant Material ..." << std::endl;

// 	voom::StVenant st(1.0, 1.0e9, 1.0/3.0);
// 	//
// 	// components of the deformation gradient tensor
// 	const double k11 = 1.6;   const double k12 = 0.0;  const double k13 = 0.0;
// 	const double k21 = 0.0;   const double k22 = 1.0;  const double k23 = 0.0;
// 	const double k31 = 0.0;   const double k32 = 0.0;  const double k33 = 1.0;	
	
// 	tvmet::Matrix<double, 3, 3>  F;
// 	F = 
// 		k11, k12, k13,
// 		k21, k22, k23,
// 		k31, k32, k33;
	
// 	st.setDeformationGradient(F);
// 	std::cout << "deformation gradient tensor = " << st.deformationGradient() << std::endl;	
// 	st.updateState(true, true, false);
// 	std::cout << "1st Piola stress = " << st.piolaStress() << std::endl;
// 	st.ConsistencyTest();
	
// 	return 0;
// }
