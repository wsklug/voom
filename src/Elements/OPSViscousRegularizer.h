/*
 * OPSViscousRegularizer.h
 *
 *  Created on: Aug 22, 2017
 *      Author: amit
 */

#ifndef SRC_ELEMENTS_OPSVISCOUSREGULARIZER_H_
#define SRC_ELEMENTS_OPSVISCOUSREGULARIZER_H_

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"
#include "Node.h"

namespace voom {

//! Simple viscous regularization
/*! ViscousRegularizer computes an energy quadratic in the difference
 between the current vector of DOF of a set of nodes and those at
 some previous state/increment/iteration (termed the reference
 state).
 This difference can be thought of as a "velocity".  The
 force is thus proportional to this velocit, hence the label
 "viscous".

 \f[
 E = \frac{k}{2}\sum_{ia} (x_{ia} - \bar{x}_{ia})^2
 \f]

 \f[
 f_{ia} = k(x_{ia}-\bar{x}_{ia})
 \f]

 \f[
 k_{iajb} = k\delta_{ab}\delta_{ij}
 \f]
 */

typedef std::vector<OPSNode*> OPSNodeContainer;
typedef std::vector<OPSNode*>::const_iterator OPSIterator;

class OPSViscousRegularizer: public Element {

public:

	//! Construct from a set of nodes and a "viscosity" proportionality factor
	OPSViscousRegularizer(const OPSNodeContainer & nodes, double viscosity) {
		//_baseNodes = nodes;
		_nodes = nodes;
		_viscosity = viscosity;
		_energy = 0.0;
		int dof = 0;
		for (OPSIterator n = _nodes.begin(); n != _nodes.end(); n++) {
			dof += 3; // OPSNode has 3 position DOFs per node
		}
		_reference.resize(dof);
		_reference = 0.0;
		step();
	}

	//! Assign the current state to the reference state
	void step() {
		int I = 0;
		for (OPSIterator n = _nodes.begin(); n != _nodes.end(); n++) {
			for (int i = 0; i < 3; i++, I++) {
				_reference(I) = (*n)->getPoint(i);
			}
		}
	}

	//! Compute the energy, force, and stiffness
	void compute(bool f0, bool f1, bool f2) {
		if (f0)
			_energy = 0.0;
		int I = 0;
		for (OPSIterator n = _nodes.begin(); n != _nodes.end();	n++) {
			for (int i = 0; i < 3; i++, I++) {
				double dx = (*n)->getPoint(i) - _reference(I);
				if (f0)
                    _energy += 0.5 * _viscosity * dx * dx;
				if (f1)
                    (*n)->addForce(i, _viscosity * dx);
			}
		}
	}

	//! return the viscosity proportionality factor
	double viscosity() const {
		return _viscosity;
	}

	//! return the viscosity proportionality factor
	void setViscosity(double v) {
		_viscosity = v;
	}

	//! compute norm of difference between reference and current
	double velocity() const {
		double v = 0.0;
		int I = 0;
		for (OPSIterator n = _nodes.begin(); n != _nodes.end(); n++) {
			for (int i = 0; i < 3; i++, I++) {
				v = std::max(v, std::abs(_reference(I) - (*n)->getPoint(i)));
			}
		}
		return v;
	}
private:

	//! proportionality factor \f$k\f$ in the energy
	double _viscosity;
	std::vector<OPSNode*> _nodes;

	//! reference state \f$\bar{x}_{ia}\f$
	blitz::Array<double, 1> _reference;
};

} // end namespace



#endif /* SRC_ELEMENTS_OPSVISCOUSREGULARIZER_H_ */
