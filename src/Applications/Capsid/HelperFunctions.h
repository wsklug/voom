#ifndef _HELPERFUNCTIONS_H
#define _HELPERFUNCTIONS_H


#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <unistd.h>
#include <tvmet/Vector.h>
#include "Node.h"
#include <limits>

#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
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

#if defined(_OPENMP)
#include <omp.h>
#endif

#if VTK_MAJOR_VERSION < 6
#define SetInputData SetInput
#endif

namespace voom
{

// Function declarations
void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
			double avgEdgeLen, double percentStrain);
void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
			double avgEdgeLen, std::vector<double> percentStrain);
void insertValenceInVtk(std::vector<std::string> fileNames);

std::vector<double> calcEdgeLenAndStdDev(std::vector<DeformationNode<3>*> a,
	std::vector< tvmet::Vector<int,3> > b);

std::vector<tvmet::Vector<int, 3> > delaunay3DSurf(std::vector<DeformationNode<3>*> a);

std::vector<std::vector<double> > getSphCellLimits(vtkSmartPointer<vtkPolyData> bins, int long_res);

void putParticlesInBins(std::vector<std::vector<double> > cellLimits,
	Eigen::Matrix3Xd newCurr, std::vector<DeformationNode<3>*> defNodes,
	vtkSmartPointer<vtkDoubleArray> binDensity, int viterMax);

/*
This function calculates radius of a LoopShellBody shell using
quadrature points
*/
template<class T>
std::vector<double> getRadialStats(LoopShellBody<T>* bd, 
	tvmet::Vector<double, 3> Xavg) {

	typedef typename LoopShellBody<T> LSB;
	typedef typename LoopShell<T> LS;

	LSB::FeElementContainer elements = bd->shells();
	std::vector<double> qpRadius(elements.size(), 0.0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int e = 0; e < elements.size(); e++) {

		const LS::NodeContainer eleNodes = elements[e]->nodes();
		LS::QuadPointContainer quadPoints = elements[e]->quadraturePoints();

		for (LS::ConstQuadPointIterator quadPoint = quadPoints.begin();
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
