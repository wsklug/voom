#ifndef _HELPERFUNCTIONS_H
#define _HELPERFUNCTIONS_H


#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <tvmet/Vector.h>
#include "Node.h"
#include <limits>
#include <set>
#include <list>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkExtractEdges.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkNew.h>
#include <vtkLine.h>
#include <vtkIdList.h>
#include <vtkUnsignedIntArray.h>
#include <vtkDataArray.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkMeshQuality.h>
#include <Eigen/Dense>
#include <LoopShellBody.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCleanPolyData.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

namespace voom
{
	//Structure declaration
	struct neighbors {
		vtkIdType _id;
		double _distance;
		neighbors(vtkIdType id, double dist) :_id(id), _distance(dist) {}
		bool operator<(const neighbors& e)const { return _distance < e._distance; }
	};

	//Unique cells
	struct uniqueTriangle {
		vtkIdType _id1, _id2, _id3;
		uniqueTriangle(vtkIdType id1, vtkIdType id2, 
			vtkIdType id3):_id1(id1), _id2(id2), _id3(id3) {}
		bool operator<(const uniqueTriangle& e)const { 
			std::vector<vtkIdType> set1 = {_id1, _id2, _id3};
			std::vector<vtkIdType> set2 = {e._id1, e._id2, e._id3};
			std::sort(set1.begin(), set1.end());
			std::sort(set2.begin(), set2.end());
			if (set1[0] != set2[0]) 
				return (set1[0] < set2[0]);
			else if (set1[1] != set2[1]) 
				return (set1[1] < set2[1]);
			else 
				return (set1[2] < set2[2]);
		}
	};

// Function declarations
void writeEdgeStrainVtk(std::vector<std::string> &fileNames, \
			double avgEdgeLen, double percentStrain);
void writeEdgeStrainVtk(std::vector<std::string> &fileNames, \
			double avgEdgeLen, std::vector<double> percentStrain);
void insertValenceInVtk(std::vector<std::string> &fileNames);

std::vector<double> calcEdgeLenAndStdDev(const std::vector<DeformationNode<3>*> &a,
	const std::vector< tvmet::Vector<int,3> > &b);

std::vector<tvmet::Vector<int, 3> > delaunay3DSurf(const std::vector<DeformationNode<3>*> &a);

void meshSphericalPointCloud(const vtkSmartPointer<vtkPolyData> pd, double searchRad, 
	const std::string fileName);

//std::vector<tvmet::Vector<int, 3> > Poisson3DSurf(const std::vector<DeformationNode<3>*> &a);

std::vector<std::vector<double> > getSphCellLimits(const vtkSmartPointer<vtkPolyData> &bins, int long_res);

void putParticlesInBins(const std::vector<std::vector<double> > &cellLimits,
	const Eigen::Matrix3Xd &newCurr, const std::vector<DeformationNode<3>*> &defNodes,
	vtkSmartPointer<vtkDoubleArray> binDensity, int viterMax);

std::vector<double> getMeshQualityInfo(const std::vector<DeformationNode<3>*> &a, 
	const std::vector< tvmet::Vector<int, 3> > connectivities);

void plotMorseBonds(const std::vector<std::string> &fileNames, std::string fname, 
	double epsilon, double Rshift, double sigma, vtkSmartPointer<vtkCellArray> bonds);

std::vector<double> getRadialStats(vtkSmartPointer<vtkPolyData> pd,
	tvmet::Vector<double, 3> Xavg);

/*
This function calculates radius of a LoopShellBody shell using
quadrature points
*/
template<class T>
std::vector<double> getRadialStats(LoopShellBody<T>* bd, 
	tvmet::Vector<double, 3> Xavg) {

	typedef LoopShellBody<T> LSB;
	typedef LoopShell<T> LS;

	typename LSB::FeElementContainer elements = bd->shells();
	std::vector<double> qpRadius(elements.size(), 0.0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int e = 0; e < elements.size(); e++) {

		const typename LS::NodeContainer eleNodes = elements[e]->nodes();
		typename LS::QuadPointContainer quadPoints = elements[e]->quadraturePoints();

		for (typename LS::ConstQuadPointIterator quadPoint = quadPoints.begin();
			quadPoint != quadPoints.end(); ++quadPoint) {

			LoopShellShape s = (*quadPoint).shape;
			const LoopShellShape::FunctionArray fn = s.functions();
			tvmet::Vector<double, 3> Xq(0.0);

			for (int i = 0; i < fn.size(); i++) {
				Xq += tvmet::mul(eleNodes[i]->point(), fn(i));
			}

			double qpR = tvmet::norm2(Xq - Xavg);
			qpRadius[e] = qpR;
		}
	}

	double Ravg = 0.0;
	for (int i = 0; i < qpRadius.size(); i++) {
		Ravg += qpRadius[i];
	}
	Ravg /= qpRadius.size();

	double dRavg2 = 0.0;
	for (int i = 0; i < qpRadius.size(); i++) {
		double dR = qpRadius[i] - Ravg;
		dRavg2 += dR*dR;
	}
	dRavg2 /= qpRadius.size();

	double asphericity = dRavg2 / (Ravg*Ravg);

	std::vector<double> output = { Ravg, asphericity };
	return output;
}

}
#endif
