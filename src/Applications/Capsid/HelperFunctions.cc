#include "HelperFunctions.h"

///////////////////////////////////////////////////////////////////////////////
//                                                                           //    
//                        WRITEEDGESTRAINVTK BEGINS                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//The method writeEdgeStrainVtk() inserts edge strain information in a vtk file
//The calling method must ensure that 'fileName' exists and is a valid vtk file.
//MUST use the extension '.vtk' in 'fileName'.
namespace voom
{
	using namespace std;

	void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
		double avgEdgeLen, double percentStrain) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < fileNames.size(); i++) {
			//for(std::vector<std::string>::iterator it=fileNames.begin();
			//    it!= fileNames.end(); ++it){
			//  std::string fileName = *it;
			std::string fileName = fileNames[i];
			//Check that the file exists
			assert(ifstream(fileName.c_str()));

			vtkSmartPointer<vtkPolyDataReader> reader =
				vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(fileName.c_str());

			reader->ReadAllScalarsOn();
			reader->ReadAllVectorsOn();

			vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
			reader->Update();

			vtkSmartPointer<vtkPolyDataWriter> writer =
				vtkSmartPointer<vtkPolyDataWriter>::New();

			//Following few lines of code are meant to obtain number of edges
			//from the mesh
			vtkSmartPointer<vtkExtractEdges> extractEdges =
				vtkSmartPointer<vtkExtractEdges>::New();
			extractEdges->SetInputConnection(reader->GetOutputPort());

			//extractEdges->Update();
			vtkSmartPointer<vtkPolyData> wireFrame = extractEdges->GetOutput();
			extractEdges->Update();

			vtkSmartPointer<vtkCellArray> lines = wireFrame->GetLines();
			int numLines = lines->GetNumberOfCells();

			//string vectorName="displacements";
			vtkSmartPointer<vtkDataArray> displacements = wireFrame->GetPointData()->
				GetVectors("displacements");

			//The following vtkDoubleArray will be used to store the
			//strain in each edge
			vtkSmartPointer<vtkDoubleArray> edgeStrain =
				vtkSmartPointer<vtkDoubleArray>::New();
			edgeStrain->SetNumberOfComponents(1);
			edgeStrain->SetNumberOfTuples(numLines);
			edgeStrain->SetName("EdgeStrains");

			vtkIdType npts;
			vtkIdType *pts;
			lines->InitTraversal();
			vtkIdType index = 0;
			while (lines->GetNextCell(npts, pts)) {
				double p1[3], p2[3];
				double disp1[3], disp2[3];
				wireFrame->GetPoint(pts[0], p1);
				wireFrame->GetPoint(pts[1], p2);
				displacements->GetTuple(pts[0], disp1);
				displacements->GetTuple(pts[1], disp2);
				tvmet::Vector<double, 3> p1v(p1[0] + disp1[0], p1[1] + disp1[1], p1[2] + disp1[2]);
				tvmet::Vector<double, 3> p2v(p2[0] + disp2[0], p2[1] + disp2[1], p2[2] + disp2[2]);
				tvmet::Vector<double, 3> line(p1v - p2v);
				double strain = (tvmet::norm2(line) - avgEdgeLen) / avgEdgeLen;
				edgeStrain->SetTuple1(index++, strain);
			}
			wireFrame->GetCellData()->SetScalars(edgeStrain);

			int numPoints = wireFrame->GetNumberOfPoints();

			//The following vtkDoubleArray will be used to store the
			//number of unstable bonds for each particle
			vtkSmartPointer<vtkUnsignedIntArray> unstableBonds =
				vtkSmartPointer<vtkUnsignedIntArray>::New();
			unstableBonds->SetNumberOfComponents(1);
			unstableBonds->SetNumberOfTuples(numPoints);
			unstableBonds->SetName("unstableBonds");

			//cellIds will be used to temporarily hold the CELLS that use a
			//point specified by a point id.
			vtkSmartPointer<vtkIdList> cellIds =
				vtkSmartPointer<vtkIdList>::New();

			wireFrame->BuildLinks();

			unsigned int weakBonds;
			double currEdgeStrain;
			vtkIdType currEdge;

			for (int p = 0; p < numPoints; p++) {
				weakBonds = 0;
				wireFrame->GetPointCells(p, cellIds);
				for (int z = 0; z < cellIds->GetNumberOfIds(); z++) {
					currEdge = cellIds->GetId(z);
					currEdgeStrain = edgeStrain->GetTuple1(currEdge);
					if (currEdgeStrain > percentStrain*0.01) {
						weakBonds++;
					}
				}
				unstableBonds->SetValue(p, weakBonds);
				cellIds->Reset();
			}

			wireFrame->GetPointData()->AddArray(unstableBonds);

			//The following array will store approximate strain in each
			//triangle
			vtkSmartPointer<vtkDoubleArray> avgStrain =
				vtkSmartPointer<vtkDoubleArray>::New();
			avgStrain->SetNumberOfComponents(1);
			avgStrain->SetNumberOfTuples(mesh->GetNumberOfCells());
			avgStrain->SetName("ApproxEleStrain");

			tvmet::Vector<double, 3> x1, x2, x3; //vertex position vectors
			double e31, e32, e12; //edge lengths

			int ntri = mesh->GetNumberOfCells();
#ifdef _OPENMP
#pragma omp parallel for private(x1,x2,x3,e31,e32,e12)
#endif    
			for (int j = 0; j < ntri; j++) {
				assert(mesh->GetCell(j)->GetNumberOfPoints() == 3);
				mesh->GetPoint(mesh->GetCell(j)->GetPointId(0), &(x1[0]));
				mesh->GetPoint(mesh->GetCell(j)->GetPointId(1), &(x2[0]));
				mesh->GetPoint(mesh->GetCell(j)->GetPointId(2), &(x3[0]));
				e31 = tvmet::norm2(x3 - x1);
				e32 = tvmet::norm2(x3 - x2);
				e12 = tvmet::norm2(x1 - x2);
				double temp = (std::abs(e31 - avgEdgeLen) +
					std::abs(e32 - avgEdgeLen) +
					std::abs(e12 - avgEdgeLen)) / (3.0*avgEdgeLen);
				avgStrain->SetValue(j, temp);
			}

			mesh->GetCellData()->AddArray(avgStrain);

			std::stringstream sstm;
			std::string tempFile;
			sstm << fileName << "-bak.vtk";
			tempFile = sstm.str();
			writer->SetFileName(tempFile.c_str());
			writer->SetInputData(mesh);
			writer->Update();
			//writer->SetFileTypeToBinary();
			writer->Write();
			if (std::remove(fileName.c_str()) != 0) {
				perror("Error deleting file:");
			}
			if (std::rename(tempFile.c_str(), fileName.c_str()) != 0) {
				perror("Error renaming file:");
			}

			/*
		  The next few lines of code involve string manipulations to comeup
		  with good file-names that Paraview can recognize to be sequence of
		  the same simulation. A file name like "T7-relaxed-10.vtk" is
		  converted to "T7-EdgeStrain-10.vtk"
			*/
			int pos = fileName.find(".vtk");
			if (pos != -1) {
				fileName.erase(pos, string::npos);
			}
			pos = fileName.find("relaxed-");
			if (pos != -1) {
				fileName.erase(pos, 8);
				pos = fileName.find("-");
				string serialNum = fileName.substr(pos + 1, string::npos);
				fileName.erase(pos, string::npos);
				fileName = "./" + fileName + "-EdgeStrain-" + serialNum + ".vtk";
			}
			else {
				fileName = "./" + fileName + "-EdgeStrain.vtk";
			}

			writer->SetFileName(fileName.c_str());
			writer->SetInputData(wireFrame);
			writer->Update();
			//writer->SetFileTypeToBinary();
			writer->Write();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//                                                                      //
	//                        OVERLOADED VERSION                            //
	//                                                                      //
	//////////////////////////////////////////////////////////////////////////

	void writeEdgeStrainVtk(std::vector<std::string> fileNames, \
		double avgEdgeLen, std::vector<double> percentStrain) {

		assert(fileNames.size() == percentStrain.size());
		for (int i = 0; i < fileNames.size(); i++) {
			std::vector<std::string> fakeVector;
			fakeVector.push_back(fileNames[i]);

			writeEdgeStrainVtk(fakeVector, avgEdgeLen, percentStrain[i]);
		}
	}

	///////////////////////////////////////////////////////////////////////////
	//                                                                       //
	//                    CALCEDGELENANDSTDDEV BEGINS                        //
	//                                                                       //
	///////////////////////////////////////////////////////////////////////////
	/*
	  Calculates average edge lengths of triangles in the mesh and the
	  standard deviation in the edge lengths.
	*/

	std::vector<double> calcEdgeLenAndStdDev
	(std::vector< DeformationNode<3>* > defNodes,
		vector< tvmet::Vector<int, 3> > connectivities) {

		double EdgeLength = 0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < connectivities.size(); i++) {
			std::vector<int> cm(3);
			for (int j = 0; j < 3; j++) cm[j] = connectivities[i](j);
			// Edge vectors in current config.
			tvmet::Vector<double, 3>
				e31(defNodes[cm[0]]->point() - defNodes[cm[2]]->point()),
				e32(defNodes[cm[1]]->point() - defNodes[cm[2]]->point()),
				e12(defNodes[cm[1]]->point() - defNodes[cm[0]]->point()),
				eCent(defNodes[cm[2]]->point());
			// Compute average edge length for each triangle
			double temp =
				(tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12)) / 3.0;

#pragma omp atomic
			EdgeLength += temp;
		}
		EdgeLength /= connectivities.size();

		// Calculate the standard deviation of side lengths of the
		// equilateral triangles
		double stdDevEdgeLen = 0.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < connectivities.size(); i++) {
			std::vector<int> cm(3);
			for (int j = 0; j < 3; j++) cm[j] = connectivities[i](j);
			// Edge vectors in current config.
			tvmet::Vector<double, 3>
				e31(defNodes[cm[0]]->point() - defNodes[cm[2]]->point()),
				e32(defNodes[cm[1]]->point() - defNodes[cm[2]]->point()),
				e12(defNodes[cm[1]]->point() - defNodes[cm[0]]->point()),
				eCent(defNodes[cm[2]]->point());
			double temp = std::pow(tvmet::norm2(e31) - EdgeLength, 2.0) +
				std::pow(tvmet::norm2(e32) - EdgeLength, 2.0) +
				std::pow(tvmet::norm2(e12) - EdgeLength, 2.0);

#pragma omp atomic
			stdDevEdgeLen += temp;
		}

		stdDevEdgeLen /= connectivities.size();
		stdDevEdgeLen = sqrt(stdDevEdgeLen);

		std::vector<double> result;
		result.push_back(EdgeLength);
		result.push_back(stdDevEdgeLen);
		return result;
	}

	///////////////////////// INSERTVALENCEINVTK BEGINS ///////////////////////////
	/////////////////////////                           //////////////////////////

	//The method insertValenceInVtk() inserts valence information in a vtk file
	//The calling method must ensure that 'fileName' exists and is a valid vtk file.
	//MUST use the extension '.vtk' in 'fileName'.
	//'mesh' is a pointer to a vtkPolyData

	void insertValenceInVtk(std::vector<std::string> fileNames) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < fileNames.size(); i++) {

			//Check that the file exists
			assert(ifstream(fileNames[i].c_str()));

			//vtkDataSetReader * reader = vtkDataSetReader::New();

			vtkSmartPointer<vtkPolyDataReader> reader =
				vtkSmartPointer<vtkPolyDataReader>::New();
			reader->SetFileName(fileNames[i].c_str());
			reader->ReadAllScalarsOn();
			reader->ReadAllVectorsOn();

			vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
			reader->Update();

			int numPoints = pd->GetNumberOfPoints();

			//The following vtkUnsignedIntArray will be used to store the
			//number of CELLS in the mesh that share the POINT denoted by the
			//index of the vector
			vtkSmartPointer<vtkUnsignedIntArray> countPointCells =
				vtkSmartPointer<vtkUnsignedIntArray>::New();
			countPointCells->SetNumberOfComponents(1);
			countPointCells->SetNumberOfTuples(numPoints);
			countPointCells->SetName("Valence");

			//cellIds will be used to temporarily hold the CELLS that use a
			//point specified by a point id.
			vtkSmartPointer<vtkIdList> cellIds =
				vtkSmartPointer<vtkIdList>::New();

			pd->BuildLinks();
			for (int p = 0; p < pd->GetNumberOfPoints(); p++) {
				pd->GetPointCells(p, cellIds);
				countPointCells->SetValue(p, cellIds->GetNumberOfIds());
				cellIds->Reset();
			}

			pd->GetPointData()->AddArray(countPointCells);
			//We will use a vtkPolyDataWriter to write our modified output
			//files that will have Capsomer information as well
			vtkSmartPointer<vtkPolyDataWriter> writer =
				vtkSmartPointer<vtkPolyDataWriter>::New();
			std::stringstream sstm;
			std::string tempFileName;
			sstm << "Temp_" << i << ".vtk";
			tempFileName = sstm.str();
			writer->SetFileName(tempFileName.c_str());
			writer->SetInputData(pd);
			//writer->SetFileTypeToBinary();
			writer->Write();
			if (std::remove(fileNames[i].c_str()) != 0) {
				perror("Error deleting file:");
			}
			if (std::rename(tempFileName.c_str(), fileNames[i].c_str()) != 0) {
				perror("Error renaming file:");
			}
		}
	}

}
