// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__PotentialBody_h__)
#define __PotentialBody_h__

#include <blitz/array.h>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "Body.h"
#include "Potential.h"
#include "PotentialElement.h"
#include "voom.h"
#include "Node.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

using namespace std;

namespace voom
{

	/*!
	  PotentialBody is a concrete class derived from Body, implementing
	  potential for a cloud of 3D points
	*/

	class PotentialBody : public Body
	{
	public:
		//! Construct body from material, nodes, volumes, and LME parameters
		PotentialBody(Potential * Mat, const vector<DeformationNode<3> * > & DefNodes,
			double SearchR);


		//! Destructor
		~PotentialBody() {
			for (uint i = 0; i < _elementVector.size(); i++) {
				delete _elementVector[i];
			}
		};

		//! Do mechanics on Body
		void compute(bool f0, bool f1, bool f2);

		void recomputeNeighbors(double searchR);

		//! Return the energy of the body
		double totalStrainEnergy() const { return _energy; };

		//! General printing of a Paraview file
		void printParaview(const string name) const;

		virtual void pushBack(Element* e) {
			_elements.push_back(e);
		}

		const vector<PotentialElement* > & getPotentialElements() {
			return _elementVector;
		}



	private:
		// Potential material
		Potential * _mat;

		// Potential elements
		vector<PotentialElement* > _elementVector;

		// Nodes
		const vector<DeformationNode<3> * > & _defNodes;

		// Search radius
		double _searchR;

#ifdef WITH_MPI
		int _processorRank;
		int _nProcessors;
#endif

	}; // Potential body

} // namespace voom

#endif // __PotentialBody_h__
