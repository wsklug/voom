// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
//

#include "RigidHemisphereAL.h"

namespace voom
{

	RigidHemisphereAL::RigidHemisphereAL(const DefNodeContainer & nodes,
		double k, double R, Vector3D xc,
		double friction) {
		_defNodes = nodes;
		for (ConstDefNodeIterator n = _defNodes.begin(); n != _defNodes.end(); n++)
			_baseNodes.push_back(*n);

		_k = k;
		_R = R;
		_xc = xc;
		_mu = friction;
		_FZ = 0.0;
		_penetration = 0.0;
		_active.resize(nodes.size());
		_forces.resize(_defNodes.size());
		_forces = Vector3D(0.0);
		_pressureMultipliers.resize(_defNodes.size());
		_pressureMultipliers = 0.0;
		_frictionMultipliers.resize(_defNodes.size());
		_frictionMultipliers = Vector3D(0.0);
		_normals.resize(_defNodes.size());
		_normals = Vector3D(0.0);
		updateContact();

	}

	int RigidHemisphereAL::active() const {
		int n = 0;
		for (int i = 0; i < _active.size(); i++) {
			if (_active[i]) n++;
		}
		return n;
	}

	void RigidHemisphereAL::updateContact() {

		// Now update surface normals and check for active contacts
		for (int a = 0; a < _defNodes.size(); a++) {
			const Vector3D & x = _defNodes[a]->point();
			double R = tvmet::norm2(x - _xc);

			// recompute outward pointing normal vector 
			Vector3D & n = _normals(a);
			n = x - _xc;
			n /= R;

			// check for active contacts
			if (R <= _R) {  // we are in contact (penetrating)	
				_active[a] = true;
			}
			else {
				_active[a] = false;
			}

			// Uzawa update: set multipliers equal to most recently computed
			// forces
			_pressureMultipliers(a) = dot(_forces(a), _normals(a));
			_frictionMultipliers(a) = _forces(a)
				- _pressureMultipliers(a)*_normals(a);


		}

		return;
	}

	void RigidHemisphereAL::compute(bool f0, bool f1, bool f2) {
		//std::cout << "RigidHemisphere::compute()" << std::endl;
		if (f0) _energy = 0.0;

		if (f1) {
			_FZ = 0.0;
			_forces = Vector3D(0.0);
		}

		_penetration = 0.0;

		for (int a = 0; a < _defNodes.size(); a++) {

			DefNode * nd = _defNodes[a];
			Vector3D x;
			x = nd->point();

			// current normal
			Vector3D n;
			n = x - _xc;
			double R = tvmet::norm2(n);
			n /= R;

			// Do normal contact first

			// normal gap function
			double gn = R - _R;
			_penetration = std::max(_penetration, -gn);

			// pressure multiplier and force
			double pn = _pressureMultipliers(a);
			double fn = pn + _k*gn;
			if (fn > 0) { // case 1: no contact 
				if (f0) {
					_energy -= 0.5*sqr(pn) / _k;
				}

			}
			else { // case 2: contact

				if (f0) {
					_energy += (pn + 0.5*_k*gn)*gn;
				}
				if (f1) {
					_forces(a) += fn*n;
				}
			}

			// now frictional contact 

			// previous contact point = xc + radius times previous normal

			// compute vector from previous contact point to current point
			Vector3D dx;
			dx = x - (_xc + _R*_normals(a));

			// tangential gap, (I-n^t n)*dx
			Vector3D gt;
			gt = dx - dot(_normals(a), dx)*_normals(a);

			// friction multiplier and force
			const Vector3D pt = _frictionMultipliers(a);
			double norm_pt = norm2(pt);

			Vector3D ft;
			ft = pt + _k*gt;

			double norm_ft = norm2(ft);

			if (pn > 0) { // case 1: no contact 
				if (f0) {
					_energy -= 0.5*sqr(norm_pt) / _k;
				}

			}
			else if (norm_ft <= _mu*std::abs(pn)) { // case 2: contact, stick

				if (f0) {
					_energy += dot(pt + 0.5*_k*gt, gt);
				}
				if (f1) {
					_forces(a) += ft - dot(ft, _normals(a))*_normals(a);
				}

			}
			else { // case 3: contact, slip

				if (f0) {
					_energy -=
						0.5*(sqr(norm_pt) - 2.0*_mu*std::abs(pn)*norm_ft + sqr(_mu*pn)) / _k;
				}
				if (f1) {
					_forces(a) +=
						_mu*std::abs(pn)*(ft - dot(ft, _normals(a))*_normals(a)) / norm_ft;
				}

			}

			if (f1) {
				const Vector3D & f = _forces(a);
				_FZ += f(2);
				for (int i = 0; i < 3; i++) nd->addForce(i, f(i));
			}

		}
		return;
	}


} // namespace voom
