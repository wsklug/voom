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
#include <blitz/array-impl.h>
#include "VoomMath.h"

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkWarpVector.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
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

		//! General printing of a Paraview file
		void printParaview(const string name, 
			Eigen::Matrix3Xd, 
			vector< tvmet::Vector<int, 3> > connectivities) const;

		virtual void pushBack(Element* e) {
			_elements.push_back(e);
		}

		const vector<PotentialElement* > & getPotentialElements() {
			return _elementVector;
		}

		//! Return the mean squared displacement
		double rmsd();

		//!Get nearest neighbor for rmsd calculation
		std::vector<int> initialNearestNeighbor();

	protected:
		// Potential material
		Potential * _mat;

		// Potential elements
		vector<PotentialElement* > _elementVector;

		// Nodes
		const vector<DeformationNode<3> * > & _defNodes;

		// Search radius
		double _searchR;

		//Nearest neighbor in Reference config, only for rmsd
		std::vector<int> _nearestNeighbor;

#ifdef WITH_MPI
		int _processorRank;
		int _nProcessors;
#endif

	}; // Potential body

} // namespace voom

#endif // __PotentialBody_h__
