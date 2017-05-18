#include "MorsePeriodic.h"

namespace voom {

	std::vector<double> MorsePeriodic::getPeriodicBox() {
		return _box;
	}

	Vector3D MorsePeriodic::getPeriodicDiffVector(DeformationNode<3>* nodeA,
		DeformationNode<3>* nodeB) {
		Vector3D diff;
		diff = nodeA->point() - nodeB->point();
		for (int k = 0; k < 3; k++) {
			if (diff(k) > _box[k] * 0.5) {
				diff(k) = diff(k) - _box[k];
			}
			if (diff(k) <= -1 * _box[k] * 0.5) {
				diff(k) = diff(k) + _box[k];
			}
		}
		return diff;
	}

	void voom::MorsePeriodic::updateState(DeformationNode<3>* nodeA,
		DeformationNode<3>* nodeB, bool fl0, bool fl1, bool fl2)
	{
		Vector3D diff;
		diff = getPeriodicDiffVector(nodeA, nodeB);
		double r = tvmet::norm2(diff);
		if (fl0) {
			// Morse energy
			_W = _epsilon*(exp(-2 * _sigma*(r - _Rshift)) - 2 * exp(-_sigma*(r - _Rshift)));
		}

		if (fl1) {
			// Morse forces
			double factor = ((2.0*_sigma*exp(-(r - _Rshift)*_sigma) - 2.0*_sigma*exp(-2.0*(r - _Rshift)*_sigma))
				*_epsilon) / r;
			Vector3D ForceIncrement(0.0);
			ForceIncrement = factor*(diff);
			nodeA->updateForce(ForceIncrement);
			ForceIncrement = -ForceIncrement;
			nodeB->updateForce(ForceIncrement);
		}

		if (fl2) {
			// Not implemented 
			cout << "Stiffness calculations in Morse potential not implemented " << endl;
			// exit(1);
		}
	}

	double voom::MorsePeriodic::computeTension(DeformationNode<3>* nodeA, DeformationNode<3>* nodeB)
	{
		Vector3D diff;
		diff = getPeriodicDiffVector(nodeA, nodeB);
		double r = tvmet::norm2(diff);
		return (2.0*_sigma*exp(-(r - _Rshift)*_sigma) - 2.0*_sigma*exp(-2.0*(r - _Rshift)*_sigma))
			*_epsilon;
	}
} // namespace voom