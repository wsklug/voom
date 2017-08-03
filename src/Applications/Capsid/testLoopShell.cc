#include <vector>
#include <string>
#include <tvmet/Vector.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataNormals.h>
#include "voom.h"
#include "Node.h"
#include "FVK.h"
#include "LoopShellBody.h"
#include "TriangleQuadrature.h"
#include "HelperFunctions.h"

using namespace tvmet;
using namespace std;
using namespace voom;

typedef FVK MaterialType;
typedef LoopShellBody<MaterialType> LSB;

int main(int argc, char* argv[] ){
    
    string inFile;
    inFile = (argc > 1)? argv[1] : "RegularPatch.vtk";    
    
    vtkSmartPointer<vtkPolyDataReader> reader =
		vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(inFile.c_str());	
    reader->Update();
    /*
    vtkSmartPointer<vtkPolyDataNormals> normals =
    vtkSmartPointer<vtkPolyDataNormals>::New();    
    normals->SetInputConnection(reader->GetOutputPort());
    normals->ComputeCellNormalsOn();
    normals->ConsistencyOn();
    normals->SplittingOff();
    normals->Update();
	vtkSmartPointer<vtkPolyData> mesh = normals->GetOutput();
	*/
    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
	std::cout << "Mesh->GetNumberOfPoints() = " << mesh->GetNumberOfPoints()
		<< std::endl;
	/*
	// create vector of nodes
	int dof = 0;
	std::vector< NodeBase* > nodes;
    
	// read in points
	for (int a = 0; a < mesh->GetNumberOfPoints(); a++) {
		int id = a;
		DeformationNode<3>::Point x;
		mesh->GetPoint(a, &(x[0]));
		NodeBase::DofIndexMap idx(3);
		for (int j = 0; j < 3; j++) idx[j] = dof++;
		DeformationNode<3>* n = new DeformationNode<3>(id, idx, x);
		nodes.push_back(n);
	}
	assert(nodes.size() != 0);

	// read in triangle connectivities
	vector< tvmet::Vector<int, 3> > connectivities;
	tvmet::Vector<int, 3> c;
	int ntri = mesh->GetNumberOfCells();
	connectivities.reserve(ntri);
	std::cout << "Number of triangles: " << ntri << endl;

	for (int i = 0; i < ntri; i++) {
		assert(mesh->GetCell(i)->GetNumberOfPoints() == 3);
		for (int a = 0; a < 3; a++) c[a] = mesh->GetCell(i)->GetPointId(a);
		connectivities.push_back(c);
	}
    
    double KC = 10.0;
    double KG = 10.0;
    double C0 = 10.0;
    
    int quadOrder = 1;
    
    MaterialType bending(KC, KG, C0, 0.0, 0.0);
    LSB * bd = new LSB(bending, connectivities, nodes, quadOrder);
    bd->setOutput(paraview);
    
    double cleanTol = 0.01;
    int loopSurfSubDiv = 0;
    vtkSmartPointer<vtkPolyData> lssPd = 
        bd->getLoopShellSurfPoints(cleanTol, loopSurfSubDiv);
    */
    meshSphericalPointCloud(mesh, 0.144, "LoopShellSurface.vtk");

    vtkSmartPointer<vtkPolyDataWriter> writer = 
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(mesh);
    writer->SetFileName("Refined.vtk");
    writer->Write();
    return(1);
}
