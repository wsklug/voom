// -*- C++ -*-
//----------------------------------------------------------------------
//
//                  William S. Klug & Luigi Perotti
//                University of California Los Angeles
//                  (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#include "PotentialBody.h"

namespace voom
{
	PotentialBody::PotentialBody(Potential * Mat, const vector<DeformationNode<3> * > & DefNodes,
		double SearchR) :
		_mat(Mat), _defNodes(DefNodes), _searchR(SearchR)
	{

#ifdef WITH_MPI
		MPI_Comm_size(MPI_COMM_WORLD, &_nProcessors);
		MPI_Comm_rank(MPI_COMM_WORLD, &_processorRank);
#endif
		_dof = 0;
		// Initialize Body.h containers
		_nodes.insert(_nodes.begin(), DefNodes.begin(), DefNodes.end());
		for (ConstNodeIterator n = _nodes.begin(); n != _nodes.end(); n++) {
			_dof += (*n)->dof();
		}

		// Initialize material objects
		_elementVector.reserve(_defNodes.size());

		for (uint i = 0; i < _defNodes.size(); i++)
		{
			Vector3D CenterNode = _defNodes[i]->point();
			set<DeformationNode<3> *> domain;
			// Find neighbors of CenterNode
			for (uint j = 0; j < _defNodes.size(); j++)
			{
				if (i != j && tvmet::norm2(CenterNode - _defNodes[j]->point()) <= _searchR) {
					domain.insert(_defNodes[j]);
				}
			}
			// Build potential element
			PotentialElement * el = new PotentialElement(_mat, _defNodes[i], domain);
			_elementVector.push_back(el);
		}

		//Initialize initial nearest neighbor vector;
		std::vector<int> temp(_defNodes.size(),-1);
		_nearestNeighbor = temp;

	}; // PotentialBody constructor

	void PotentialBody::recomputeNeighbors(const double searchR) {
		_searchR = searchR;

		for (uint i = 0; i < _defNodes.size(); i++)
		{
			Vector3D CenterNode = _defNodes[i]->point();
			set<DeformationNode<3> *> domain;
			// Find neighbors of CenterNode
			for (uint j = 0; j < _defNodes.size(); j++)
			{
				if (i != j && tvmet::norm2(CenterNode - _defNodes[j]->point())
					<= _searchR) {
					domain.insert(_defNodes[j]);
				}
			}
			// Modify domain of each element
			_elementVector[i]->resetDomain(domain);

		}

	}	

	//! Compute E0, E1, E2
	void PotentialBody::compute(bool f0, bool f1, bool f2)
	{
		// Initialize energy to be zero
		if (f0) _energy = 0.0;

		// Loop over material objects
		for (uint i = 0; i < _elementVector.size(); i++)
		{
			_elementVector[i]->compute(f0, f1, f2);
			_energy += _elementVector[i]->energy();
		}

		//If the inherited _elements container is non-empty
		for (uint i = 0; i < _elements.size(); i++)
		{
			_elements[i]->compute(f0, f1, f2);
			_energy += _elements[i]->energy();
		}


		return;
	};

	//! General printing of a Paraview file
	void PotentialBody::printParaview(const string name) const {

		std::string fileName = name + ".vtk";
		vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> pd
			= vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points
			= vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> displacements =
			vtkSmartPointer<vtkDoubleArray>::New();

		tvmet::Vector<double, 3> refPosition, currPosition, disp;
		points->SetNumberOfPoints(_defNodes.size());
		displacements->SetNumberOfComponents(3);
		displacements->SetNumberOfTuples(_defNodes.size());
		displacements->SetName("displacements");

		for (int i = 0; i < _defNodes.size(); i++) {
			double X[3] = { 0.0,0.0,0.0 };
			refPosition = _defNodes[i]->position();
			currPosition = _defNodes[i]->point();
			disp = currPosition - refPosition;
			for (int q = 0; q < 3; q++) X[q] = refPosition(q);
			points->SetPoint(i, X);
			displacements->SetTuple3(i, disp(0), disp(1), disp(2));
		}
		pd->SetPoints(points);
		pd->GetPointData()->AddArray(displacements);

		writer->SetInputData(pd);
		writer->SetFileName(fileName.c_str());
		writer->Write();

	}

	void PotentialBody::printParaview(const string name,
		Eigen::Matrix3Xd transformed, 
		vector< tvmet::Vector<int, 3> > connectivities) const {
		std::string fileName = name + ".vtk";
		vtkSmartPointer<vtkPolyDataWriter> writer =
			vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> pd
			= vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPoints> points
			= vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkDoubleArray> displacements =
			vtkSmartPointer<vtkDoubleArray>::New();

		tvmet::Vector<double, 3> refPosition, currPosition, disp;
		points->SetNumberOfPoints(_defNodes.size());
		displacements->SetNumberOfComponents(3);
		displacements->SetNumberOfTuples(_defNodes.size());
		displacements->SetName("displacements");

		for (int i = 0; i < _defNodes.size(); i++) {
			double X[3] = { 0.0,0.0,0.0 };
			refPosition = _defNodes[i]->position();
			for (int q = 0; q < 3; q++) {
				X[q] = refPosition(q);
				currPosition(q) = transformed(q,i);
			}
			disp = currPosition - refPosition;
			points->SetPoint(i, X);			
			displacements->SetTuple3(i, 
				disp(0), disp(1), disp(2));
		}
		
		//Prepare the Cell data
		vtkSmartPointer<vtkCellArray> cells =
			vtkSmartPointer<vtkCellArray>::New();
		for (int i = 0; i < connectivities.size(); i++) {
			vtkIdType cell[3];
			for (int j = 0; j < 3; j++) {				
				cell[j] = connectivities[i](j);
			}
			cells->InsertNextCell(3, cell);
		}

		pd->SetPoints(points);
		pd->GetPointData()->AddArray(displacements);
		pd->SetPolys(cells);

		writer->SetInputData(pd);
		writer->SetFileName(fileName.c_str());
		writer->Write();
	}

	//! Mean Sqaured Displacement
	double PotentialBody::rmsd() {
		double rmsd = 0; //root Mean Squared Displacement
		int nn = -1;
		for (int i = 0; i < _defNodes.size(); i++) {
			Vector3D xi, xj, diff;
			nn = _nearestNeighbor[i];
			xi = (_defNodes[i]->point() - _defNodes[i]->position());
			xj = (_defNodes[nn]->point() - _defNodes[nn]->position());
			diff = xi - xj;
			rmsd += tvmet::dot(xi, xi);
		}
		rmsd = rmsd / (2 * _defNodes.size());
		return rmsd;
	}

	//! Initial nearest neighbors, used in rmsd calculations
	std::vector<int> PotentialBody::initialNearestNeighbor()
	{
		//Identify nearest neighbor
		Vector3D xi, xj;
		for (int p = 0; p < _defNodes.size(); p++) {
			xi = _defNodes[p]->position();
			double dist_min = 10e6;
			for (int q = 0; q < _defNodes.size(); q++) {
				if (p != q) {
					xj = _defNodes[q]->position();
					double dist = tvmet::norm2(xi - xj);
					if (dist < dist_min) {
						dist_min = dist;
						_nearestNeighbor[p] = q;
					}
				}
			}
		}
		return _nearestNeighbor;
	}

} // namespace voom
