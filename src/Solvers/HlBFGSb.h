// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                University of California Los Angeles
//                   (C) 2006 All Rights Reserved
//
//----------------------------------------------------------------------
//
//
//----------------------------------------------------------------------
//
// Reference:
//  Headers only implementation of lBFGSb
// https://yixuan.cos.name/LBFGSpp/
//
/////////////////////////////////////////////////////////////////////////

/*!
\file HlBFGSb.h

\brief Interface to a concrete class for L-BFGS-B solver for static
equilibrium of a (nonlinear) Finite Element model.

*/

#if !defined(__HlBFGSb_h__)
#define __HlBFGSb_h__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <blitz/array.h>
#include <vector>
#include "Lbfgsb.h"
#include <Eigen/Core>
#include <LBFGS.h>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace LBFGSpp;
using namespace boost::multiprecision;
using namespace Eigen;

namespace voom
{
	//template < class Scalar >
	class HlBFGSb
	{
	public:
		typedef Eigen::Matrix< cpp_dec_float_50, Eigen::Dynamic, 1> VectorT;
		typedef blitz::Array<double, 1> Vector_t;

		//! Default Constructor
		HlBFGSb(int n, Model& model, bool d = false) : _model(model), 
			_debug(d)
		{
			_x.resize(n,1);
			_g.resize(n,1);
		}

		//! destructor
		virtual ~HlBFGSb() {}

		void zeroOutData(bool f0, bool f1, bool f2) {
			if (f0) _f = 0.0;
			//if (f1) _g.col(0).setZero();
			if (f1) _g = 0.0;
		}

		double & field(int i) { return _x(i); }
		double & function() { return _f; }
		double & gradient(int i) { return _g(i); }
		double & hessian(int i, int j) {
			std::cerr << "No stiffness in Lbfgsb solver." << std::endl;
			exit(0);
		}
		const double hessian(int i, int j) const {
			std::cerr << "No stiffness in Lbfgsb solver." << std::endl;
			exit(0);
		};

		cpp_dec_float_50 operator()(VectorT& _x, VectorT& _g){			
			// copy starting guess from model
			_model.getField(*this);
			_computeAll();
			return _f;
		}

		double & hessian(int i) { return hessian(i, i); }
		const double hessian(int i) const { return hessian(i, i); }
		const blitz::Array<double, 2> & hessian() const;

	private:
		Model& _model;
		double _f;
		Vector_t _x;
		Vector_t _g;
		double _factr, _pgtol;
		bool _debug;

		void _computeAll() {
			_model.putField(*this);
			_model.computeAndAssemble(*this, true, true, false);
			return;
		}

	};
} // namespace voom

#endif // __HlBFGSb_h__
