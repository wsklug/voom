/*! 
 \file ContractionWrapper.h
 */

#include "Material.h"
#include "Math.h"

#if !defined(__ContractionWrapper_h__)
#define __ContractionWrapper_h__

namespace voom {

template<class Material_T>

class ContractionWrapper : public Material {
protected:
	Material_T _myMaterial;
	Tensor3D _CC; // Green deformation tensor C = trans(F)*F
	double _T;

public:

	ContractionWrapper() {}
	//ContractionWrapper(const Material_T &);
	ContractionWrapper(const Material_T &Input) {
		//if (this._myMaterial == &Input)
		//	return;
		_myMaterial = Input;
	}

	void setDeformationGradient( const Tensor3D & F ){
		std::cout << "In wrapper: F(I,J) = " << std::endl;
		std::cout << F << std::endl;
	    if( determinant(F) > 0 ) {
	      _F=F;
	      _myMaterial.setDeformationGradient(F);
	      _CC = (tvmet::trans(_F))*_F;
	    } else{
	      std::cout << "The determinant of deformation gradient tensor is not greater than zero!" << " J = " << determinant(F) << std::endl;
	    }
	  };
	void updateState(bool f0, bool f1, bool f2) {
		_myMaterial.updateState(f0, f1, f2);
		_P = _myMaterial.piolaStress();
		_C = _myMaterial.cauchyStress();
		_W = _myMaterial.energyDensity();
		std::cout << "P =" << std::endl;
		std::cout << _P << std::endl;
		// modify _P and _W 
		_P += _T*_CC*_F;
		std::cout << _P << std::endl;
		_W += (1/2)*_T*tvmet::trace(_CC*_CC);
	}

	double getT() const {
		return _T;
	}
	void setT(double t) {
		_T = t;
	}

	//void ConsistencyTest();
	void ConsistencyTest() {
		std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
		updateState(true, true, false);
		Tensor3D PAna, PNum;
		PAna = _P; // current 1st Piola stress
		const double eps = 1.0e-8 * max(_F); // perturbation
		for (int i = 0; i < PAna.rows(); i ++) {
			for (int j = 0; j < PAna.cols(); j++) {
				_F(i, j) += eps;
				setDeformationGradient(_F);
				updateState(true, true, false);
				double W = _W;
				_F(i, j) -= 2*eps;
				setDeformationGradient(_F);
				updateState(true, true, false);
				W -= _W;
				_F(i, j) += eps; // restore value
				setDeformationGradient(_F);
				PNum(i, j) = W/2.0/eps;
			}
		}
		std::cout << "ContractionWrapper.h: Analytical value of 1st Piola stress:"
				<< std::endl;
		std::cout << PAna << std::endl;
		std::cout << "ContractionWrapper.h: Numberical value of 1st Piola stress:"
				<< std::endl;
		std::cout << PNum << std::endl;
		std::cout << "ContractionWrapper.h: Difference between 1st Piola stresses:"
				<< std::endl;
		PNum -= PAna;
		PNum /= max(PAna);
		std::cout << PNum << std::endl;

		std::cout << "\n\n\n\n" << std::endl;
		return;
	}

};

//ContractionWrapper::ContractionWrapper(const Material_T &Input) {
//	if (this._myMaterial == &Input)
//		return;
//	_myMaterial = new Material_T(&Input);
//}

//void ContractionWrapper::ConsistencyTest() {
//	std::cout << "checking consistency of 1st Piola stress tensor" << std::endl;
//	updateState(true, true, false);
//	Tensor3D PAna, PNum;
//	PAna = _P; // current 1st Piola stress
//	const double eps = 1.0e-8 * max(_F); // perturbation
//	for (int i = 0; i < PAna.rows(); i ++) {
//		for (int j = 0; j < PAna.cols(); j++) {
//			_F(i, j) += eps;
//			updateState(true, true, false);
//			double W = _W;
//			_F(i, j) -= 2*eps;
//			updateState(true, true, false);
//			W -= _W;
//			_F(i, j) += eps; // restore value
//			PNum(i, j) = W/2.0/eps;
//		}
//	}
//	std::cout << "ContractionWrapper.h: Analytical value of 1st Piola stress:"
//			<< std::endl;
//	std::cout << PAna << std::endl;
//	std::cout << "ContractionWrapper.h: Numberical value of 1st Piola stress:"
//			<< std::endl;
//	std::cout << PNum << std::endl;
//	std::cout << "ContractionWrapper.h: Difference between 1st Piola stresses:"
//			<< std::endl;
//	PNum -= PAna;
//	PNum /= max(PAna);
//	std::cout << PNum << std::endl;
//
//	std::cout << "\n\n\n\n" << std::endl;
//	return;
//}
}

#endif //  !defined(__ContractionWrapper_h__)
