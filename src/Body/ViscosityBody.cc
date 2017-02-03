// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
#include <string>
#include <fstream>
#include <blitz/array-impl.h>
#include "VoomMath.h"
#include "ViscosityBody.h"


#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{
	ViscosityBody::ViscosityBody(const vector<DeformationNode<3>* > & DefNodes,
		const vector<tvmet::Vector<int, 3> > & Connectivity,
		double sigma) :
		_defNodes(DefNodes), _sigma(sigma)
	{

#ifdef WITH_MPI
		MPI_Comm_size(MPI_COMM_WORLD, &_nProcessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &_processorRank);
#endif

		// Initialize Body.h containers
		_nodes.insert(_nodes.begin(), DefNodes.begin(), DefNodes.end());
		for (ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
			_dof += (*n)->dof();
		}

		// Temporary vector to store node pairs
		vector<vector<uint > > SpringsNodesID;
		vector<uint > ind(4, 0);
		ind[1] = 1; ind[2] = 2;
		for (uint i = 0; i < Connectivity.size(); i++) {
			for (uint j = 0; j < 3; j++) {
				bool IsThere = false;
				vector<uint > spring(2, 0);
				spring[0] = Connectivity[i][ind[j]];
				spring[1] = Connectivity[i][ind[j + 1]];

				// Check if spring element already present
				for (uint k = 0; k < SpringsNodesID.size(); k++) {
					if ((SpringsNodesID[k][0] == spring[0] && SpringsNodesID[k][1] == spring[1]) ||
						(SpringsNodesID[k][0] == spring[1] && SpringsNodesID[k][1] == spring[0]))
					{
						IsThere = true;
						break;
					}
				}

				if (IsThere == false) {
					SpringsNodesID.push_back(spring);
					vector<DeformationNode<3> *> OneSpringNodes;
					OneSpringNodes.push_back(_defNodes[spring[0]]);
					OneSpringNodes.push_back(_defNodes[spring[1]]);

					SpringsNodes.push_back(OneSpringNodes);

					double OneSpringLength = tvmet::norm2(_defNodes[spring[0]]->point() - _defNodes[spring[1]]->point());
					SpringsLengths.push_back(OneSpringLength);
				}

			}
		}

	}; // Viscosity body constructor



	//! Compute E0, E1, E2
	void ViscosityBody::compute(bool f0, bool f1, bool f2)
	{
		// Energy
		if (f0) {
			_energy = 0.0;
			for (uint i = 0; i < SpringsNodes.size(); i++)
			{
				double r = tvmet::norm2((SpringsNodes[i][0])->point() - (SpringsNodes[i][1])->point());
				_energy += 0.5*_sigma*pow(r - SpringsLengths[i], 2.0);
			}
		}

		if (f1) {
			for (uint i = 0; i < SpringsNodes.size(); i++)
			{
				double r = tvmet::norm2((SpringsNodes[i][0])->point() - (SpringsNodes[i][1])->point());

				double factor = _sigma*(r - SpringsLengths[i]) / r;
				Vector3D ForceIncrement(0.0);
				ForceIncrement = factor*((SpringsNodes[i][0])->point() - (SpringsNodes[i][1])->point());

				(SpringsNodes[i][0])->updateForce(ForceIncrement);
				ForceIncrement *= -1;
				(SpringsNodes[i][1])->updateForce(ForceIncrement);
			}
		}

		return;
	};



	void ViscosityBody::resetRefConf()
	{
		for (uint i = 0; i < SpringsNodes.size(); i++)
		{
			SpringsLengths[i] = tvmet::norm2((SpringsNodes[i][0])->point() - (SpringsNodes[i][1])->point());
		}
	}


} // namespace voom
