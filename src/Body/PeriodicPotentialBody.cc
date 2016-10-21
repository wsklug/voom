// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------
#include "PeriodicPotentialBody.h"

namespace voom {

	//! Assuming a rectangular brick shaped bounding box for Periodic Boundary
	// Lx-by-Ly-by-Lz centered at the origin
	void PeriodicPotentialBody::recomputeNeighbors(double searchR)
	{
		double Lx = _boundingBox[0];
		double Ly = _boundingBox[1];
		double Lz = _boundingBox[2];

		if (Lx*0.5 < searchR || Ly*0.5 < searchR || Lz*0.5 < searchR) {
			std::cout << "Bounding box should have dimensions greater"
				<< "than twice the search radius." << std::endl;
			exit(EXIT_FAILURE);
		}
		_searchR = searchR;
		for (uint i = 0; i < _defNodes.size(); i++)
		{
			Vector3D CenterNode = _defNodes[i]->point();
			set<DeformationNode<3> *> domain;
			double distance = 0;
			Vector3D dX;
			// Find neighbors of CenterNode
			for (uint j = 0; j < _defNodes.size(); j++)
			{
				dX = CenterNode - _defNodes[j]->point();
				for (int k = 0; k < 3; k++) {
					if (dX(k) > _boundingBox[k] * 0.5) {
						dX(k) = dX(k) - _boundingBox[k];
					}
					if (dX(k) <= -1 * _boundingBox[k] * 0.5) {
						dX(k) = dX(k) + _boundingBox[k];
					}
				}
				distance = tvmet::norm2(dX);
				if (i != j && distance <= _searchR) {
					domain.insert(_defNodes[j]);
				}
			}
			// Modify domain of each element
			_elementVector[i]->resetDomain(domain);
		}
	}//Periodic Boundary Neighbors

	 //! Re-enter particles that cross a periodic boundary
	/*
	if (rx[i]<0.0) { rx[i]+=L; ix[i]--; }
	if (rx[i]>L)   { rx[i]-=L; ix[i]++; }
	*/
	void PeriodicPotentialBody::adjustPositions() {
		std::vector<double> box = _boundingBox;

		for (int i = 0; i < _defNodes.size(); i++) {
			Vector3D pos = _defNodes[i]->point();
			for (int k = 0; k < 3; k++) {
				if (pos(k) < -0.5*box[k]){
					pos(k) = pos(k) + box[k];
					_boundaryCrossCounter[i][k]--;
				}
				if (pos(k) > 0.5*box[k]) {
					pos(k) = pos(k) - box[k];
					_boundaryCrossCounter[i][k]++;
				}
			}
		}

	}

	//! Mean Sqaured Displacement
	double PeriodicPotentialBody::meanSquaredDisplacement() {
		double msd = 0; //Mean Squared Displacement
		for (int i = 0; i < _defNodes.size(); i++) {
			Vector3D tempDisp;
			tempDisp = (_defNodes[i]->point() - _defNodes[i]->position());
			for (int k = 0; k < 3; k++) {
				tempDisp(k) = tempDisp(k) + 
					_boundaryCrossCounter[i][k] * _boundingBox[k];
			}
			msd += tvmet::dot(tempDisp, tempDisp);
		}
	}

} // namespace voom