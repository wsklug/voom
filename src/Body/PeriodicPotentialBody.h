// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__PeriodicPotentialBody_h__)
#define __PeriodicPotentialBody_h__

#include "PotentialBody.h"
#include <vtkIntArray.h>

using namespace std;

namespace voom {
	class PeriodicPotentialBody : public PotentialBody {
	public:
		//! Construct body from material, nodes, volumes, and LME parameters
		PeriodicPotentialBody(Potential * Mat, const vector<DeformationNode<3> * > & DefNodes,
			double SearchR, std::vector<double> boundingBox);

		//! Destructor
		~PeriodicPotentialBody() {};

		//! Recompute neighbors for Periodic Boundary conditions
		// as per the minimum image convention
		void recomputeNeighbors(double searchR);

		//! Re-enter any particle that moves out of the box.
		void adjustPositions();

		//! Return the mean squared displacement
		double rmsd();

		//! General printing of a Paraview file
		void printParaview(const string name) const;

		//Accessor for bounding box
		std::vector<double> getBoundingBox() {
			return _boundingBox;
		}

		//Accessor for boundary cross counter
		std::vector<std::vector<int> > getBoundaryCrossCount() {
			return _boundaryCrossCounter;
		}

	private:
		std::vector<double> _boundingBox;
		std::vector<std::vector<int> > _boundaryCrossCounter;
		std::vector<int> _nearestNeighbor;
	};
} // namespace voom
#endif // __PeriodicPotentialBody_h__