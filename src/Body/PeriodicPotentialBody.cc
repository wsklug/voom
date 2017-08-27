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
	PeriodicPotentialBody::PeriodicPotentialBody(Potential * Mat, 
		const vector<DeformationNode<3>*>& DefNodes, double SearchR,
		std::vector<double> boundingBox): PotentialBody(Mat, DefNodes, SearchR)
	{
		_boundingBox = boundingBox;
		std::vector<vector<int> > A(DefNodes.size(), vector<int>(3));
		_boundaryCrossCounter = A;		
	}

	//! Assuming a rectangular brick shaped bounding box for Periodic Boundary
	// Lx-by-Ly-by-Lz centered at the origin
	void PeriodicPotentialBody::recomputeNeighbors(double searchR)
	{
		double Lx = _boundingBox[0];
		double Ly = _boundingBox[1];
		double Lz = _boundingBox[2];

		if (Lx*0.5 < searchR || Ly*0.5 < searchR || Lz*0.5 < searchR) {
			std::cout << "Bounding box should have dimensions greater"
				<< " than twice the search radius." << std::endl;
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
		int crossedCount = 0;
		bool crossedFlag = false;
		std::vector<double> box = _boundingBox;
		for (int i = 0; i < _defNodes.size(); i++) {
			crossedFlag = false;
			Vector3D pos = _defNodes[i]->point();
			for (int k = 0; k < 3; k++) {
				if (pos(k) <= -0.5*box[k]){
					pos(k) = pos(k) + box[k];
					_boundaryCrossCounter[i][k]--;
					crossedFlag = true;
				}
				if (pos(k) > 0.5*box[k]) {
					pos(k) = pos(k) - box[k];
					_boundaryCrossCounter[i][k]++;
					crossedFlag = true;
				}
			}
			crossedCount = crossedFlag ? crossedCount + 1 : crossedCount;
			if(crossedFlag) _defNodes[i]->setPoint(pos);
		}
		if (crossedCount > 0)
			std::cout << crossedCount << " particles crossed boundary."
			<< std::endl;
	}

	//! Mean Sqaured Displacement
	double PeriodicPotentialBody::rmsd() {
		double rmsd = 0; //root Mean Squared Displacement
		int nn = -1;
		for (int i = 0; i < _defNodes.size(); i++) {
			Vector3D xi,xj,diff;
			nn = _nearestNeighbor[i];
			xi = (_defNodes[i]->point() - _defNodes[i]->position());
			xj = (_defNodes[nn]->point() - _defNodes[nn]->position());
			for (int k = 0; k < 3; k++) {
				xi(k) = xi(k) +
					_boundaryCrossCounter[i][k] * _boundingBox[k];
				xj(k) = xj(k) +
					_boundaryCrossCounter[nn][k] * _boundingBox[k];
			}
			diff = xi - xj;
			rmsd += tvmet::dot(diff, diff);
		}
		rmsd = sqrt(rmsd / (2*_defNodes.size()));
		return rmsd;
	}

	void PeriodicPotentialBody::printParaview(const string name) const {
		std::string fileName = name + ".vtk";
		vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> pd
			= vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points
			= vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> displacements =
			vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkIntArray> boundaryCrossCount =
			vtkSmartPointer<vtkIntArray>::New();

		tvmet::Vector<double, 3> refPosition, currPosition, disp;
		points->SetNumberOfPoints(_defNodes.size());
		displacements->SetNumberOfComponents(3);
		displacements->SetNumberOfTuples(_defNodes.size());
		displacements->SetName("displacements");
		boundaryCrossCount->SetNumberOfComponents(3);
		boundaryCrossCount->SetNumberOfTuples(_defNodes.size());
		boundaryCrossCount->SetName("boundaryCrossCount");

		for (int i = 0; i < _defNodes.size(); i++) {
			double X[3] = { 0.0,0.0,0.0 };
			refPosition = _defNodes[i]->position();
			currPosition = _defNodes[i]->point();
			disp = currPosition - refPosition;
			for (int q = 0; q < 3; q++) {
				X[q] = refPosition(q);
			}
			points->SetPoint(i, X);
			displacements->SetTuple3(i, disp(0), disp(1), disp(2));
			boundaryCrossCount->SetTuple3(i, 
				_boundaryCrossCounter[i][0], _boundaryCrossCounter[i][1],
				_boundaryCrossCounter[i][2]);
		}
		pd->SetPoints(points);
		pd->GetPointData()->AddArray(displacements);
		pd->GetPointData()->AddArray(boundaryCrossCount);

		writer->SetInputData(pd);
		writer->SetFileName(fileName.c_str());
		writer->Write();
	}

} // namespace voom
