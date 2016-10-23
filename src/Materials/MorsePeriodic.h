// -*- C++ -*-
//----------------------------------------------------------------------
//
//                   William S. Klug, Luigi Perotti
//                University of California Los Angeles
//                 (C) 2004-2016 All Rights Reserved
//
//----------------------------------------------------------------------
//
/*!
\file Morse.h
\brief Interface for Morse potential in 3D cartesian space for Periodic
\boundary condition
*/
#ifndef _MORSEPERIODIC_H_
#define _MORSEPERIODIC_H_

#include "Morse.h"

using namespace std;

namespace voom {

	class MorsePeriodic : public Morse
	{
	public:

		// Constructors/destructors:
		//! Default constructor
		MorsePeriodic() : Morse(0.0, 0.0, 0.0), _box(3, 0) {};
		//! Construct potential from necessary constants
		MorsePeriodic(double epsilon, double sigma, double Rshift,
			const std::vector<double> bounds) : Morse(epsilon, sigma, Rshift) {
			_box = bounds;
		}

		//! Destructor
		virtual ~MorsePeriodic() {}
		
		//! Accessor for Periodic boundary
		std::vector<double> getPeriodicBox();

		//! Mutator for Periodic boundary
		void setPeriodicBox(std::vector<double> newBox) { _box = newBox; }

		void updateState(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB, bool f0, bool f1, bool f2);
		double computeTension(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB);

	protected:
		Vector3D getPeriodicDiffVector(DeformationNode<3> *nodeA, DeformationNode<3> *nodeB);

	private:
		//! Bounds of the orthorhombic periodic boundary box
		std::vector<double> _box;

	};

}// namespace voom

#endif // _MORSEPERIODIC_H_