// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug & Feng Feng
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
 \file Node.h

 \brief Node is a class for a Finite Element node representing a
 point in a mesh where some field is interpolated.

 */

#ifndef _NODE_
#define _NODE_

#ifdef _X
#undef _X
#endif
#ifdef _A
#undef _A
#endif
#ifdef _B
#undef _B
#endif

#include <vector>
#include <fstream>
#include "NodeBase.h"

namespace voom {

//! A node with point representing a value of a scalar field at
//  specified position.
template<int dim_n>
class ScalarFieldNode: public NodeBase {

public:
	typedef NodeBase Base;
	typedef double Point;
	typedef typename tvmet::Vector<double, dim_n> PositionVector;

	//! default constructor
	ScalarFieldNode() {
	}

	//! construct from position and point
	ScalarFieldNode(int id, const NodeBase::DofIndexMap & index,
			const PositionVector & X, const Point & p) :
				Base(id, index) {
		_X = X;
		_point = p;
	}

	ScalarFieldNode(int id, const NodeBase::DofIndexMap & index,
			const PositionVector & X) :
				Base(id, index) {
		_X = _point = X;
	}

	//! access reference position
	const PositionVector & position() {
		return _X;
	}

	//! access point
	const Point & point() const {
		return _point;
	}

	//! access force
	const Point & force() const {
		return _force;
	}

	//! set force
	void setForce(double f) {
		_force = f;
	}

	//! add force
	void addForce(double df) {
#ifdef _OPENMP
#pragma omp critical
#endif
		_force = _force + df;
	}

	virtual void setPosition(int i, double x) {
		assert(i < dim_n);
		_X(i) = x;
	}

	virtual double getPosition(int i) const {
		assert(i < dim_n);
		return _X(i);
	}

	// it seems like this should be inherited from Node, but GCC objects...?
	virtual void setPoint(const Point & p) {
		_point = p;
	}

	virtual double getPoint(int i) const {
		assert(i == 0);
		return _point;
	}

	virtual void setPoint(int i, double x) {
		assert(i == 0);
		_point = x;
	}

	virtual void addPoint(int i, double dx) {
		assert(i == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
		_point = _point + dx;
	}

	virtual double getForce(int i) const {
		assert(i == 0);
		return _force;
	}

	virtual void setForce(int i, double f) {
		assert(i == 0);
		_force = f;
	}

	virtual void addForce(int i, double df) {
		assert(i == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
		_force = _force + df;
	}

	virtual void setPosition(const PositionVector & p) {
		_X = p;
	}

	//! update force by some increment
	virtual void updateForce(const Point & f) {
#ifdef _OPENMP
#pragma omp critical
#endif
		_force = _force + f;
	}

	void resetPosition() {
		_X = _point;
	}

	virtual int dof() const {
		return 1;
	}

protected:
	Point _point;
	Point _force;
	PositionVector _X;
};

//! A node with point and position of equal dimension.
/*!
 A DeformationNode is a PositionNode for which the point is a
 vector of the same dimension as the position.  This is typical for
 a finite element nodes interpolating displacements or the
 deformation mapping of a bulk solid.
 */
template<int dim_n>
class DeformationNode: public NodeBase {

public:
	typedef NodeBase Base;
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n> PositionVector;

	//! default constructor
	DeformationNode() {
	}

	//! construct from position and point
	DeformationNode(int id, const NodeBase::DofIndexMap & index,
			const Point & X, const Point & p) :
				Base(id, index) {
		_X = X;
		_point = p;
	}

	DeformationNode(int id, const NodeBase::DofIndexMap & index,
			const Point & X) :
				Base(id, index) {
		_X = _point = X;
	}

	//! access reference position
	const PositionVector & position() {
		return _X;
	}

	void setPosition(const PositionVector & p) {
		_X = p;
	}

	//! access point
	const Point & point() const {
		return _point;
	}

	//! access force
	const Point & force() const {
		return _force;
	}

	virtual void setPosition(int i, double x) {
		assert(i < dim_n);
		_X(i) = x;
	}

	virtual double getPosition(int i) const {
		assert(i < dim_n);
		return _X(i);
	}

	// it seems like this should be inherited from Node, but GCC objects...?
	virtual void setPoint(const Point & p) {
		_point = p;
	}

	virtual double getPoint(int i) const {
		assert(i < dim_n);
		return _point(i);
	}

	virtual void setPoint(int i, double x) {
		assert(i < dim_n);
		_point(i) = x;
	}

	virtual void addPoint(int i, double dx) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_point(i) = _point(i) + dx; //+= dx;
	}

	virtual double getForce(int i) const {
		assert(i < dim_n);
		return _force(i);
	}

	virtual void setForce(int i, double f) {
		assert(i < dim_n);
		_force(i) = f;
	}

	virtual void addForce(int i, double df) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_force(i) = _force(i) + df; // += df;
	}

	//! update force by some increment
	virtual void updateForce(const Point & f) {
		for (int i = 0; i < dim_n; i++) {
#ifdef _OPENMP
#pragma omp critical
#endif
			_force(i) = _force(i) + f(i);
		}
	}

	void resetPosition() {
		_X = _point;
	}

	virtual int dof() const {
		return dim_n;
	}

	virtual double getStiffness(int i) const {
		assert(i < dim_n);
		return _stiff(i);
	}

	virtual void setStiffness(int i, double k) {
		assert(i < dim_n);
		_stiff(i) = k;
	}

	virtual void addStiffness(int i, double dk) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_stiff(i) = _stiff(i) + dk; // += df;
	}

	//! update stiffness by some increment
	virtual void updateStiffness(const Point & k) {
		for (int i = 0; i < dim_n; i++) {
#ifdef _OPENMP
#pragma omp critical
#endif
			_stiff(i) = _stiff(i) + k(i);
		}
	}

protected:
	Point _point;
	Point _force;
	PositionVector _X;

	Point _stiff;
};

//! A DeformationNode with its point projected along a unit normal.
/*!  The idea of this node is to allow the point to change only
 along a specified normal direction.  To elements and bodies it
 should appear as a typical DeformationNode, with a vector point
 and force.  But to the Model (and Solver) it should only have one
 dof.  Model and solver only use the methods getPoint(i),
 setPoint(i,p), addPoint(i,dp), and getForce(i).  So each of these
 needs to act like there is only one dof and one force component;
 these are the projections along the normal direction of the vector
 point and force.  Other methods that change these vectors need not
 be modified since they are used only by elements (and maybe
 bodies?).
 */
template<int dim_n>
class ProjectedDeformationNode: public DeformationNode<dim_n> {

public:
	typedef DeformationNode<dim_n> Base;
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n> PositionVector;

	using Base::_X;
	using Base::_point;
	using Base::_force;

	//! default constructor
	ProjectedDeformationNode(int id, const NodeBase::DofIndexMap & index) :
		Base(id, index) {
		const Point zero(0.0);
		_n = zero;
	}

	//! yet another constructor
	ProjectedDeformationNode(int id, const NodeBase::DofIndexMap & index,
			const Point & X, const Point & p, const Point & n) :
				DeformationNode<dim_n>(id, index, X, p) {
		_n = n;
		if (norm2(n) > 0.0)
			_n /= norm2(n);
	}
	;

	//! destructor
	virtual ~ProjectedDeformationNode() {
	}

	//! access normal
	const Point & normal() const {
		return _n;
	}
	//! access a component of normal
	const double normal(int i) const {
		assert(i < dim_n);
		return _n(i);
	}
	//! assign normal
	virtual void setNormal(const Point & n) {
		_n = n;
		if (norm2(n) > 0.0)
			_n /= norm2(n);
	}

	//! One dof (projection of point onto normal)
	int dof() const {
		return 1;
	}

	virtual double getPoint(int i) const {
		assert(i == 0);
		return dot(_point, _n);
	}

	virtual void setPoint(int i, double x) {
		assert(i == 0);
		_point = _point - _n * dot(_point, _n) + x * _n;
	}

	virtual void addPoint(int i, double dx) {
		assert(i == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
		_point = _point + _n * dx;
	}

	virtual double getForce(int i) const {
		assert(i == 0);
		return dot(_force, _n);
	}

	virtual void setForce(int i, double f) {
		assert(i < dim_n);
		_force(i) = f;
	}

	virtual void addForce(int i, double df) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_force(i) = _force(i) + df; // += df;
	}

protected:
	Point _n;

};

//! A DeformationNode for dynamics.
/*!
 This node has in addition to point and position, vectors for
 velocity and acceleration, i.e., the first and second time derivatives
 of point.  (All four quantities are of the same dimension.)
 */
template<int dim_n>
class MotionNode: public DeformationNode<dim_n> {
public:
	typedef DeformationNode<dim_n> Base;
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n> PositionVector;

	//! default constructor
	MotionNode(int id, const NodeBase::DofIndexMap & index) :
		Base(id, index) {
		const Point zero(0.0);
		_v = _a = zero;
	}
	//! construct from position and point
	MotionNode(int id, const NodeBase::DofIndexMap & index, const Point & X,
			const Point & p) :
				DeformationNode<dim_n>(id, index, X, p) {
		_v = _a = Point(0.0);
	}
	;

	//! construct from point
	MotionNode(int id, const NodeBase::DofIndexMap & index, const Point & p) {
		MotionNode(id, index, p, p);
	}
	;

	//! yet another constructor
	MotionNode(int id, const NodeBase::DofIndexMap & index, const Point & X,
			const Point & p, const Point & v, const Point & a) :
				DeformationNode<dim_n>(id, index, X, p) {
		_v = v;
		_a = a;
	}
	;

	//! one more constructor
	MotionNode(int id, const NodeBase::DofIndexMap & index, const Point & p,
			const Point & v, const Point & a) {
		MotionNode(id, index, p, p, v, a);
	}

	//! destructor
	virtual ~MotionNode() {
	}

	//! access velocity
	const Point & velocity() const {
		return _v;
	}
	//! access a component of velocity
	const double velocity(int i) const {
		assert(i < dim_n);
		return _v(i);
	}
	//! assign velocity
	virtual void setVelocity(const Point & v) {
		_v = v;
	}

	//! access acceleration
	const Point & acceleration() const {
		return _a;
	}
	//! access a component of acceleration
	const double acceleration(int i) const {
		assert(i < dim_n);
		return _a(i);
	}
	//! assign acceleration
	virtual void setAcceleration(const Point & a) {
		_a = a;
	}

	//! access mass
	double getMass() const {
		return _M;
	}
	//! assign mass
	void setMass(double m) {
		_M = m;
	}
	//! add mass
	void addMass(double m) {
		_M += m;
	}

protected:
	Point _v;
	Point _a;
	double _M; // nodal mass
};

//! A DeformationNode for brownian dynamics.
/*!
 This node has in addition to point and position, a vector for
 velocity and a matrix for mobility.
 */
template<int dim_n>
class BrownianNode: public DeformationNode<dim_n> {
public:
	typedef DeformationNode<dim_n> Base;
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n> PositionVector;

	typedef typename tvmet::Matrix<double, dim_n, dim_n> Matrix;

	//! default constructor
	BrownianNode(int id, const NodeBase::DofIndexMap & index) :
		Base(id, index), _v(0.0), _M(0.0), _D(0.0) {
	}
	//! construct from position and point
	BrownianNode(int id, const NodeBase::DofIndexMap & index, const Point & X,
			const Point & p) :
				DeformationNode<dim_n>(id, index, X, p), _v(0.0), _M(0.0), _D(0.0) {
		_v = Point(0.0);
	}
	;

	//! yet another constructor
	BrownianNode(int id, const NodeBase::DofIndexMap & index, const Point & X,
			const Point & p, const Point & v) :
				DeformationNode<dim_n>(id, index, X, p), _v(v), _M(0.0), _D(0.0) {
	}
	;

	//! destructor
	virtual ~BrownianNode() {
	}

	//! access velocity
	const Point & velocity() const {
		return _v;
	}
	//! access a component of velocity
	const double velocity(int i) const {
		assert(i < dim_n);
		return _v(i);
	}
	//! assign velocity
	virtual void setVelocity(const Point & v) {
		_v = v;
	}

	//! access mobility
	const Matrix & mobility() const {
		return _M;
	}
	//! access a component of mobility
	const double getMobility(int i, int j) const {
		assert(i < dim_n && j < dim_n);
		return _M(i, j);
	}

	//! assign mobility
	virtual void setMobility(const Matrix & m) {
		_M = m;
	}

	//! add mobility
	virtual void addMobility(const Matrix & m) {
		_M += m;
	}

	void addMobility(int i, int j, double m) {
		assert(i < dim_n && j < dim_n);
		_M(i, j) += m;
	}

	//! access drag
	const Matrix & drag() const {
		return _D;
	}
	//! access a component of drag
	const double getDrag(int i, int j) const {
		assert(i < dim_n && j < dim_n);
		return _D(i, j);
	}

	//! assign drag
	virtual void setDrag(const Matrix & d) {
		_D = d;
	}

	//! add drag
	virtual void addDrag(const Matrix & d) {
		_D += d;
	}

	void addDrag(int i, int j, double d) {
		assert(i < dim_n && j < dim_n);
		_D(i, j) += d;
	}

protected:
	Point _v;
	Matrix _M;
	Matrix _D; //Drag matrix
};

class MultiplierNode: public NodeBase {
public:

	typedef double Point;
	//! default constructor
	MultiplierNode(int id, const NodeBase::DofIndexMap & index) :
		NodeBase(id, index), _point(0.0) {
		;
	}

	//! another constructor
	MultiplierNode(int id, const NodeBase::DofIndexMap & index, const Point & p) :
		NodeBase(id, index), _point(p) {
	}
	;

	//! destructor
	virtual ~MultiplierNode() {
	}

	double point() const {
		return _point;
	}

	double force() const {
		return _force;
	}

	int dof() const {
		return 1;
	}
	double getPoint(int i) const {
		assert(i == 0);
		return _point;
	}
	void setPoint(int i, double x) {
		assert(i == 0);
		_point = x;
	}
	void addPoint(int i, double dx) {
		assert(i == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
		_point = _point + dx; //+= dx;
	}
	double getForce(int i) const {
		assert(i == 0);
		return _force;
	}
	void setForce(int i, double f) {
		assert(i == 0);
		_force = f;
	}
	void addForce(int i, double df) {
		assert(i == 0);
#ifdef _OPENMP
#pragma omp critical
#endif
		_force = _force + df; // += df;
	}

	// it seems like this should be inherited from Node, but GCC objects...?
	void setPoint(double p) {
		_point = p;
	}

protected:
	Point _point;
	Point _force;

};

//! A node with a point which is a vector.
/*!
 A VectorNode is a node for which the point is a
 vector of the templated dimension dim_n.
 */
template<int dim_n>
class VectorNode: public NodeBase {

public:

	typedef NodeBase Base;
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n> Vector;

	//! construct from position and point
	VectorNode(int id, const NodeBase::DofIndexMap & index, const Point & p) :
		Base(id, index) {
		_point = p;
	}

	double getPoint(int i) const {
		assert(i < dim_n);
		return _point(i);
	}

	void setPoint(int i, double x) {
		assert(i < dim_n);
		_point(i) = x;
	}

	void addPoint(int i, double dx) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_point(i) = _point(i) + dx; //+= dx;
	}

	double getForce(int i) const {
		assert(i < dim_n);
		return _force(i);
	}

	void setForce(int i, double f) {
		assert(i < dim_n);
		_force(i) = f;
	}

	void addForce(int i, double df) {
		assert(i < dim_n);
#ifdef _OPENMP
#pragma omp critical
#endif
		_force(i) = _force(i) + df; // += df;
	}

	int dof() const {
		return dim_n;
	}

protected:
	Point _point;
	Point _force;
};

template<int dim_n>
class XCNode: public NodeBase {
public:
	typedef typename tvmet::Vector<double, dim_n> Point;
	typedef typename tvmet::Vector<double, dim_n + 1> PointPlusOne;

	XCNode(int id, const NodeBase::DofIndexMap& index, const Point& X,
			const double& C) :
				NodeBase(id, index) {
		//int size=index.size();
		NodeBase::DofIndexMap xIndex(dim_n), cIndex(1);

		for (int i = 0; i < dim_n; i++)
			xIndex[i] = index[i];

		cIndex[0] = index[dim_n];

		DeformationNode<dim_n>* tempXNode = new DeformationNode<dim_n>(id,
				xIndex, X);
		MultiplierNode* tempCNode = new MultiplierNode(id, cIndex, C);

		xNode = tempXNode;
		cNode = tempCNode;

	}

	XCNode(int id, const NodeBase::DofIndexMap& index, const Point& X,
			const Point& p, const double& C) :
				NodeBase(id, index) {
		//int size=index.size();
		NodeBase::DofIndexMap xIndex(dim_n), cIndex(1);

		for (int i = 0; i < dim_n; i++)
			xIndex[i] = index[i];

		cIndex[0] = index[dim_n];

		DeformationNode<dim_n>* tempXNode = new DeformationNode<dim_n>(id,
				xIndex, X, p);
		MultiplierNode* tempCNode = new MultiplierNode(id, cIndex, C);

		xNode = tempXNode;
		cNode = tempCNode;

	}

	double getPoint(int i) const {
		if (i < dim_n)
			return xNode->getPoint(i);
		else
			return cNode->getPoint(0);
	}

	void setPoint(int i, double x) {
		if (i < dim_n)
			xNode->setPoint(i, x);
		else
			cNode->setPoint(0, x);
	}

	void addPoint(int i, double dx) {
		if (i < dim_n)
			xNode->addPoint(i, dx);
		else
			cNode->addPoint(0, dx);
	}

	double getForce(int i) const {
		if (i < dim_n)
			return xNode->getForce(i);
		else
			return cNode->getForce(0);
	}

	void setForce(int i, double f) {
		if (i < dim_n)
			xNode->setForce(i, f);
		else
			cNode->setForce(0, f);
	}

	void addForce(int i, double df) {
		if (i < dim_n)
			xNode->addForce(i, df);
		else
			cNode->addForce(0, df);
	}

	//used for update viscosity force in TwoPhaseBody.cc
	void updateForce(const Point & f) {
		xNode->updateForce(f);
	}

	const PointPlusOne force() const {
		PointPlusOne tempForce;

		for (int i = 0; i < dim_n + 1; i++) {
			tempForce(i) = getForce(i);
		}
		return tempForce;
	}

	void setPoint(const Point& p) {
		xNode->setPoint(p);
	}

	void setPosition(const Point& p) {
		xNode->setPosition(p);
	}

	void resetPosition() {
		xNode->resetPosition();
	}

	int dof() const {
		return dim_n + 1;
	}

	const Point& point() const {
		return xNode->point();
	}

	const Point& position() const {
		return xNode->position();
	}

	double concentration() const {
		return getPoint(dim_n);
	}

	double reactionCoordinate() const {
		return getPoint(dim_n);
	}

protected:
	DeformationNode<dim_n>* xNode;
	MultiplierNode* cNode;

};

/**
 * Oriented Particle System node is a node with 3-displacement and
 * 3-orientation degrees of freedom. Thus, it represents a 3D point
 * with an associated normal in space.
 */
class OPSNode: public NodeBase {
public:
	typedef tvmet::Vector<double, 6> Vector6D;

	typedef tvmet::Vector<double, 4> Quaternion;

	static Vector3D conjugation(Quaternion q, Vector3D p) {
		// We will assume that q and p are quaternions
		double q0,q1,q2,q3 ;
		double p1,p2,p3;
		q0 = q[0]; q1 = q[1]; q2 = q[2]; q3 = q[3];
		p1 = p[0]; p2 = p[1]; p3 = p[2];
		Vector3D ans(0.0);
		ans[0] = p1*q0*q0 + p1*q1*q1 - p1*q2*q2 - p1*q3*q3 - 2*p2*q0*q3 +
				2*p2*q1*q2 + 2*p3*q0*q2 + 2*p3*q1*q3;
		ans[1] = 2*p1*q0*q3 + 2*p1*q1*q2 + p2*q0*q0 - p2*q1*q1 + p2*q2*q2 -
				p2*q3*q3 - 2*p3*q0*q1 + 2*p3*q2*q3;
		ans[2] = -2*p1*q0*q2 + 2*p1*q1*q3 + 2*p2*q0*q1 + 2*p2*q2*q3 +
				p3*q0*q0 - p3*q1*q1 - p3*q2*q2 + p3*q3*q3;
		return ans;
	}

	//! Functions to convert from point normal to rotation vector or vice-versa
	static Vector3D convertRotVecToNormal(Vector3D r) {
		// Assume z-axis of Global Coord Sys is the reference for ptNormal rotation
		Vector3D zaxis(0.0);
		zaxis[2] = 1.0;
		Vector3D ptNormal(0.0);
		double vi = tvmet::norm2(r); // Angle of rotation
		// Check 0-angle case
		if (vi < 1e-10) {
			ptNormal[2] = 1.0;
		}
		else if (std::abs(vi - M_PI) < 1e-10) {
			ptNormal[2] = -1.0;
		}
		else {
			// Make unit quaternion with angle vi and axis along rotVec
			Vector3D w;
			w = (sin(0.5 * vi) / vi) * r;
			Quaternion q(cos(0.5 * vi), w[0], w[1], w[2]);
			ptNormal = conjugation(q, zaxis);
		}
		return ptNormal;
	}
	static Vector3D convertNormalToRotVec(Vector3D n) {
		double angle;
		Vector3D axis(0.0);
		// Cross-product of n with z-axis.
		Vector3D cross_prod(-n[1], n[0], 0.0 );
		double cross_prod_norm = tvmet::norm2( cross_prod );
		// Check if n is parallel or anti-parallel to z-axis
		double p3 = n[2];
		if( cross_prod_norm < 1e-10){
			axis = 1.0, 0.0, 0.0; // Arbitrarily choose the x-axis
			angle = (p3 > 0.0)? 0.0 : M_PI;
		}
		else{
			angle = asin( cross_prod_norm );
			angle = (p3 < 0.0)? (M_PI - angle) : angle;
			axis = cross_prod/cross_prod_norm;
		}
		Vector3D rotVec;
		rotVec = angle*axis;
		return rotVec;
	}

	//! construct from given reference and deformed positions
	OPSNode(int id, const NodeBase::DofIndexMap & index, const Vector3D &X,
			const Vector3D &x) :
				NodeBase(id, index) {
		Vector3D N,n,R,r;
		N = X / tvmet::norm2(X);
		R = convertNormalToRotVec(N);
		n = x / tvmet::norm2(x);
		r = convertNormalToRotVec(n);

		_X = X[0], X[1], X[2], R[0], R[1], R[2];
		_x = x[0], x[1], x[2], r[0], r[1], r[2];
	}

	OPSNode(int id, const NodeBase::DofIndexMap & index, const Vector3D & X):
		NodeBase(id, index) {
		Vector3D N, R;
		N = X / tvmet::norm2(X);
		R = convertNormalToRotVec(N);
		_X = X[0], X[1], X[2], R[0], R[1], R[2];
		_x = _X;
	}

	OPSNode(int id, const NodeBase::DofIndexMap & index, const Vector3D &X,
			const Vector3D &x, const Vector3D &N, const Vector3D &n):
				NodeBase(id, index){
		Vector3D R, r;
		R = convertNormalToRotVec(N);
		r = convertNormalToRotVec(n);
		_X = X[0], X[1], X[2], R[0], R[1], R[2];
		_x = x[0], x[1], x[2], r[0], r[1], r[2];
	}

	//! access reference position
	const Vector6D & referenceDOFs() {
		return _X;
	}

	void setReferenceDOFs(const Vector6D & p) {
		_X = p;
	}

	//! access deformed position
	const Vector6D & deformedDOFs() {
		return _x;
	}

	void setDeformedDOFs(const Vector6D & p) {
			_x = p;
	}

	//! access force
	const Vector6D & force() const {
		return _force;
	}

	void setReferenceDOF(int i, double x) {
		assert(i < 6);
		_X(i) = x;
	}

	double getReferenceDOF(int i) const {
		assert(i < 6);
		return _X(i);
	}

	double getDeformedDOF(int i) const {
		assert(i < 6);
		return _x(i);
	}

	void setDeformedDOF(int i, double x) {
		assert(i < 6);
		_x(i) = x;
	}

	void addToDeformedDOF(int i, double dx) {
		assert(i < 6);
		_x(i) = _x(i) + dx;
	}

	virtual double getPoint(int i) const{
		return getDeformedDOF(i);
	}
	virtual void setPoint(int i, double x){
		setDeformedDOF(i,x);
	}
	virtual void addPoint(int i, double dx){
		addToDeformedDOF(i,dx);
	}

	double getForce(int i) const {
		assert(i < 6);
		return _force(i);
	}

	void setForce(int i, double f) {
		assert(i < 6);
		_force(i) = f;
	}

	void addForce(int i, double df) {
		assert(i < 6);
		_force(i) = _force(i) + df;
	}

	Vector3D deformedPosition(){
		Vector3D x(_x[0],_x[1],_x[2]);
		return x;
	}

	Vector3D referencePosition(){
		Vector3D x(_X[0],_X[1],_X[2]);
		return x;
	}

	Vector3D deformedRotationVector(){
		Vector3D x(_x[3],_x[4],_x[5]);
		return x;
	}

	Vector3D referenceRotationVector(){
			Vector3D x(_X[3],_X[4],_X[5]);
			return x;
	}

	void setReferencePosAndRotVec(Vector3D X){
		Vector3D N, R;
		N = X / tvmet::norm2(X);
		R = convertNormalToRotVec(N);
		for(int i=0; i < 3; i++) _X(i) = X(i);
		for(int i=0; i < 3; i++) _X(i+3) = R(i);
	}

	void setDeformedPosAndRotVec(Vector3D X){
		Vector3D N, R;
		N = X / tvmet::norm2(X);
		R = convertNormalToRotVec(N);
		for(int i=0; i < 3; i++) _x(i) = X(i);
		for(int i=0; i < 3; i++) _x(i+3) = R(i);
	}

	//! update force by some increment
	void updateForce(const Vector6D & f) {
		for (int i = 0; i < 6; i++) {
			_force(i) = _force(i) + f(i);
		}
	}

	void resetDeformedToReference() {
		_X = _x;
	}

	int dof() const {
		return 6;
	}

protected:
	// first three DOFs are position coordinates and next 3
	//are rotation vector
	Vector6D _x;
	Vector6D _X;
	Vector6D _force;
};

}
#endif // _NODE_
