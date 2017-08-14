#include "HelperFunctions.h"

/*
 The method writeEdgeStrainVtk() inserts edge strain information in a vtk file
 The calling method must ensure that 'fileName' exists and is a valid vtk file.
 MUST use the extension '.vtk' in 'fileName'.
 */
namespace voom {
using namespace std;

void writeEdgeStrainVtk(std::vector<std::string> &fileNames, double avgEdgeLen,
		double percentStrain) {

#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (int i = 0; i < fileNames.size(); i++) {
		std::string fileName = fileNames[i];
		//Check that the file exists
		assert(ifstream(fileName.c_str()));

		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
				vtkPolyDataReader>::New();
		reader->SetFileName(fileName.c_str());

		reader->ReadAllScalarsOn();
		reader->ReadAllVectorsOn();

		vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
		reader->Update();

		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<
				vtkPolyDataWriter>::New();

		//Following few lines of code are meant to obtain number of edges
		//from the mesh
		vtkSmartPointer<vtkExtractEdges> extractEdges = vtkSmartPointer<
				vtkExtractEdges>::New();
		extractEdges->SetInputConnection(reader->GetOutputPort());

		//extractEdges->Update();
		vtkSmartPointer<vtkPolyData> wireFrame = extractEdges->GetOutput();
		extractEdges->Update();

		vtkSmartPointer<vtkCellArray> lines = wireFrame->GetLines();
		int numLines = lines->GetNumberOfCells();

		//string vectorName="displacements";
		vtkSmartPointer<vtkDataArray> displacements =
				wireFrame->GetPointData()->GetVectors("displacements");

		//The following vtkDoubleArray will be used to store the
		//strain in each edge
		vtkSmartPointer<vtkDoubleArray> edgeStrain = vtkSmartPointer<
				vtkDoubleArray>::New();
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
			tvmet::Vector<double, 3> p1v(p1[0] + disp1[0], p1[1] + disp1[1],
					p1[2] + disp1[2]);
			tvmet::Vector<double, 3> p2v(p2[0] + disp2[0], p2[1] + disp2[1],
					p2[2] + disp2[2]);
			tvmet::Vector<double, 3> line(p1v - p2v);
			double strain = (tvmet::norm2(line) - avgEdgeLen) / avgEdgeLen;
			edgeStrain->SetTuple1(index++, strain);
		}
		wireFrame->GetCellData()->SetScalars(edgeStrain);

		int numPoints = wireFrame->GetNumberOfPoints();

		//The following vtkDoubleArray will be used to store the
		//number of unstable bonds for each particle
		vtkSmartPointer<vtkUnsignedIntArray> unstableBonds = vtkSmartPointer<
				vtkUnsignedIntArray>::New();
		unstableBonds->SetNumberOfComponents(1);
		unstableBonds->SetNumberOfTuples(numPoints);
		unstableBonds->SetName("unstableBonds");

		//cellIds will be used to temporarily hold the CELLS that use a
		//point specified by a point id.
		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

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
				if (currEdgeStrain > percentStrain * 0.01) {
					weakBonds++;
				}
			}
			unstableBonds->SetValue(p, weakBonds);
			cellIds->Reset();
		}

		wireFrame->GetPointData()->AddArray(unstableBonds);

		//The following array will store approximate strain in each
		//triangle
		vtkSmartPointer<vtkDoubleArray> avgStrain = vtkSmartPointer<
				vtkDoubleArray>::New();
		avgStrain->SetNumberOfComponents(1);
		avgStrain->SetNumberOfTuples(mesh->GetNumberOfCells());
		avgStrain->SetName("ApproxEleStrain");

		tvmet::Vector<double, 3> x1, x2, x3;           //vertex position vectors
		double e31, e32, e12;                               //edge lengths

		int ntri = mesh->GetNumberOfCells();
#ifdef _OPENMP
# pragma omp parallel for private(x1,x2,x3,e31,e32,e12)
#endif
		for (int j = 0; j < ntri; j++) {
			assert(mesh->GetCell(j)->GetNumberOfPoints() == 3);
			mesh->GetPoint(mesh->GetCell(j)->GetPointId(0), &(x1[0]));
			mesh->GetPoint(mesh->GetCell(j)->GetPointId(1), &(x2[0]));
			mesh->GetPoint(mesh->GetCell(j)->GetPointId(2), &(x3[0]));
			e31 = tvmet::norm2(x3 - x1);
			e32 = tvmet::norm2(x3 - x2);
			e12 = tvmet::norm2(x1 - x2);
			double temp = (std::abs(e31 - avgEdgeLen)
					+ std::abs(e32 - avgEdgeLen) + std::abs(e12 - avgEdgeLen))
					/ (3.0 * avgEdgeLen);
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
		} else {
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

void writeEdgeStrainVtk(std::vector<std::string> &fileNames, double avgEdgeLen,
		std::vector<double> percentStrain) {

	assert(fileNames.size() == percentStrain.size());
	for (int i = 0; i < fileNames.size(); i++) {
		std::vector<std::string> fakeVector;
		fakeVector.push_back(fileNames[i]);

		writeEdgeStrainVtk(fakeVector, avgEdgeLen, percentStrain[i]);
	}
}

/*
 Calculates average edge lengths of triangles in the mesh and the
 standard deviation in the edge lengths.
 */

std::vector<double> calcEdgeLenAndStdDev(
		const std::vector<DeformationNode<3>*> &defNodes,
		const vector<tvmet::Vector<int, 3> > &connectivities) {

	double EdgeLength = 0.0;
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (int i = 0; i < connectivities.size(); i++) {
		std::vector<int> cm(3);
		for (int j = 0; j < 3; j++)
			cm[j] = connectivities[i](j);
		// Edge vectors in current config.
		tvmet::Vector<double, 3> e31(
				defNodes[cm[0]]->point() - defNodes[cm[2]]->point()), e32(
				defNodes[cm[1]]->point() - defNodes[cm[2]]->point()), e12(
				defNodes[cm[1]]->point() - defNodes[cm[0]]->point()), eCent(
				defNodes[cm[2]]->point());
		// Compute average edge length for each triangle
		double temp =
				(tvmet::norm2(e31) + tvmet::norm2(e32) + tvmet::norm2(e12))
						/ 3.0;

#pragma omp atomic
		EdgeLength += temp;
	}
	EdgeLength /= connectivities.size();

	// Calculate the standard deviation of side lengths of the
	// equilateral triangles
	double stdDevEdgeLen = 0.0;
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (int i = 0; i < connectivities.size(); i++) {
		std::vector<int> cm(3);
		for (int j = 0; j < 3; j++)
			cm[j] = connectivities[i](j);
		// Edge vectors in current config.
		tvmet::Vector<double, 3> e31(
				defNodes[cm[0]]->point() - defNodes[cm[2]]->point()), e32(
				defNodes[cm[1]]->point() - defNodes[cm[2]]->point()), e12(
				defNodes[cm[1]]->point() - defNodes[cm[0]]->point()), eCent(
				defNodes[cm[2]]->point());
		double temp = std::pow(tvmet::norm2(e31) - EdgeLength, 2.0)
				+ std::pow(tvmet::norm2(e32) - EdgeLength, 2.0)
				+ std::pow(tvmet::norm2(e12) - EdgeLength, 2.0);

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
/*
 This function returns does Delaunay3D triangulation of a cloud of points which
 form a closed shell by first projecting it to a unit sphere. We extract the
 surface polydata connectivity from the resulting shell. We return the connectivity
 after mapping them back to the original point ids
 */
vector<tvmet::Vector<int, 3> > delaunay3DSurf(
		const std::vector<DeformationNode<3>*> &nodes) {
	clock_t t1, t2;
	t1 = clock();
	//Convert DeformationNodes to vtkPoints projected to a unit sphere
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < nodes.size(); i++) {
		tvmet::Vector<double, 3> cp = nodes[i]->point();
		cp /= norm2(cp);
		pts->InsertNextPoint(cp(0), cp(1), cp(2));
	}
	//We need to insert a point at origin for 3D Delaunay triangulation
	pts->InsertNextPoint(0.0, 0.0, 0.0);
	pd->SetPoints(pts);
	//Now we will do the Delaunay triangulation
	vtkSmartPointer<vtkDelaunay3D> d3D = vtkSmartPointer<vtkDelaunay3D>::New();
	d3D->SetInputData(pd);
	d3D->Update();
	//Extract the surface from convex hull obtained after Delaunay3D
	vtkSmartPointer<vtkDataSetSurfaceFilter> dss = vtkSmartPointer<
			vtkDataSetSurfaceFilter>::New();
	dss->SetInputConnection(d3D->GetOutputPort());
	vtkSmartPointer<vtkPolyData> surf = dss->GetOutput();
	dss->Update();

	//Sanity check: Number of points in 'surf' should be equal to number
	//of points in 'surf'
	if (surf->GetNumberOfPoints() != nodes.size()) {
		exit(EXIT_FAILURE);
	}
	//We need to map the point ids of the extracted surface 'surf'
	//to the point ids of 'pd'. Both are renormalized to the surface
	//of a unit sphere. The point ids of 'pd' are same as vector index
	//of 'nodes'
	std::map<int, int> nodeIndexMap;
	for (int i = 0; i < surf->GetNumberOfPoints(); i++) {
		tvmet::Vector<double, 3> y(0.0);
		surf->GetPoint(i, &y[0]);
		for (int j = 0; j < pd->GetNumberOfPoints(); j++) {
			tvmet::Vector<double, 3> currNode(0.0);
			pd->GetPoint(j, &currNode[0]);
			if (tvmet::norm2(y - currNode) < 1e-6) {
				nodeIndexMap[i] = j;
				break;
			}
		}
	}
	//Map the connectivity information from point ids
	//of 'surf' to 'node' indices
	tvmet::Vector<int, 3> c;
	int ntri = surf->GetNumberOfCells();
	vector<tvmet::Vector<int, 3> > connectivities;
	connectivities.reserve(ntri);
	for (int i = 0; i < ntri; i++) {
		for (int a = 0; a < 3; a++) {
			c[a] = nodeIndexMap[surf->GetCell(i)->GetPointId(a)];
		}
		connectivities.push_back(c);
	}
	t2 = clock();
	float diff = ((float) t2 - (float) t1);
	std::cout << "\tRe-meshing processing time =" << diff / CLOCKS_PER_SEC
			<< " seconds" << std::endl;
	return connectivities;
}

/*
 Overloaded version of delaunay3DSurf(). It prints surface to a file and uses a
 vtkPolyData as an input.
 */
void meshSphericalPointCloud(const vtkSmartPointer<vtkPolyData> input,
		double searchRad, const std::string fileName) {
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<
			vtkPolyDataWriter>::New();
	int numPoints = input->GetNumberOfPoints();
	// Construct K-d tree to search for neighbors
	vtkSmartPointer<vtkKdTree> kdt = vtkSmartPointer<vtkKdTree>::New();
	vtkSmartPointer<vtkPoints> pts = input->GetPoints();
	kdt->BuildLocatorFromPoints(pts);

	std::vector<std::vector<vtkIdType> > capsomers;
	for (int i = 0; i < numPoints; i++) {
		tvmet::Vector<double, 3> centerNode(0.0);
		input->GetPoint(i, &centerNode[0]);
		vtkSmartPointer<vtkIdList> neighbors =
				vtkSmartPointer<vtkIdList>::New();
		kdt->FindPointsWithinRadius(searchRad, &centerNode[0], neighbors);
		// Separate centerNode from neighbors
		std::vector<vtkIdType> currCapso;
		currCapso.push_back(i);
		//std::cout<<"Neighbors of point "<<i<<std::endl;
		for (int j = 0; j < neighbors->GetNumberOfIds(); j++) {
			vtkIdType currId = neighbors->GetId(j);
			if (currId != i) {
				//std::cout<<"\t\t"<<currId<<std::endl;
				currCapso.push_back(currId);
			}
		}
		capsomers.push_back(currCapso);
	}

	/*
	 * Now we have all neighbors. We need to convert them into triangles.
	 * Form a vector with centerNode and first Neighbor. Then calculate
	 * angles of vectors formed by centernode and subsequent neighbors with
	 * the first vector using atan2. Sort by angles.
	 */
	std::set<uniqueTriangle> uniqueCells;
	for (int i = 0; i < capsomers.size(); i++) {
		std::vector<vtkIdType> currCapso = capsomers[i];
		std::set<neighbors> angles;
		tvmet::Vector<double, 3> pt0(0.0), pt1(0.0), vec0;
		if (currCapso.size() < 7) {
			std::cout << "Pentamer (or less) detected! Point id = " << i
					<< std::endl;
			std::cout << "\tNeighbors = " << std::endl;
			for (int z = 1; z < currCapso.size(); z++) {
				std::cout << currCapso[z] << " " << std::endl;
			}
		}

		vtkIdType vert0 = currCapso[0];
		vtkIdType vert1 = currCapso[1];	//TODO: Guard against seg fault
		input->GetPoint(vert0, &pt0[0]);
		input->GetPoint(vert1, &pt1[0]);
		vec0 = pt1 - pt0;
		double vec0_norm;
		vec0_norm = tvmet::norm2(vec0);
		vec0 = vec0 / vec0_norm;
		//std::cout<<" vec0 : "<< vec0[0] <<","<<vec0[1]<<","<<vec0[2] << std::endl;
		for (int j = 2; j < currCapso.size(); j++) {
			vtkIdType currVert = currCapso[j];
			tvmet::Vector<double, 3> ptj(0.0), currVec, currCross, axis;
			double currNorm, currAngle, currSin, currCos, sign;
			input->GetPoint(currVert, &ptj[0]);
			currVec = ptj - pt0;
			currNorm = tvmet::norm2(currVec);
			currVec = currVec / currNorm;
			//std::cout<<"\tvec"<< j <<" = "<<currVec[0] <<","<<currVec[1]<<","<<currVec[2] << std::endl;
			currCross = tvmet::cross(vec0, currVec);
			currSin = tvmet::norm2(currCross);
			axis = currCross / currSin;
			//std::cout<<"\taxis"<< j <<" = "<<axis[0] <<","<<axis[1]<<","<<axis[2] << std::endl;
			sign = tvmet::dot(axis, pt0);
			currSin = (sign > 0.0) ? currSin : -1.0 * currSin;
			currCos = tvmet::dot(vec0, currVec);
			//std::cout<<"\tsin =" << currSin << std::endl;
			//std::cout<<"\tcos =" << currCos << std::endl;
			currAngle = (180 / M_PI) * atan2(currSin, currCos);
			currAngle = (currAngle < 0) ? (360 + currAngle) : currAngle;
			//std::cout<<"\tangle =" << currAngle << std::endl;
			neighbors currNeighbor(currVert, currAngle);
			angles.insert(currNeighbor);
		}
		//Now we have angular positions of all neighbors. So make triangles.
		vtkIdType vert_next = vert1;
		std::set<neighbors>::iterator it;
		for (it = angles.begin(); it != angles.end(); ++it) {
			neighbors n = *it;
			vtkIdType currVert = n._id;
			uniqueTriangle cell(vert0, vert_next, currVert);
			//std::cout<<"\t\tCell = "<< vert0 <<","<<vert_next<<","<<currVert<<std::endl;
			uniqueCells.insert(cell);
			vert_next = currVert;
		}
		uniqueTriangle cell(vert0, vert_next, vert1);
		uniqueCells.insert(cell);
	}
//Write the uniqueCells to a file
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	std::set<uniqueTriangle>::iterator it;
	for (it = uniqueCells.begin(); it != uniqueCells.end(); ++it) {
		uniqueTriangle currTri = *it;
		vtkIdType cell[3] = { currTri._id1, currTri._id2, currTri._id3 };
		cells->InsertNextCell(3, cell);
	}
	input->SetPolys(cells);
//Clear the vertices from 'input'
	vtkSmartPointer<vtkCellArray> emptyVertices =
			vtkSmartPointer<vtkCellArray>::New();
	input->SetVerts(emptyVertices);
//Finally, print the 'input' with triangles
	writer->SetFileName(fileName.c_str());
	writer->SetInputData(input);
	writer->Write();
}

/*
 This method creates a 3D surface from unorganized points using
 Poisson surface reconstruction. It returns a connectivity matrix
 */
/*
 vector<tvmet::Vector<int, 3> > Poisson3DSurf(const std::vector<DeformationNode<3>*>
 &nodes)
 {
 vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
 for (int i = 0; i < nodes.size(); i++) {
 tvmet::Vector<double, 3> cp = nodes[i]->point();
 pts->InsertNextPoint(cp(0), cp(1), cp(2));
 }
 vtkSmartPointer<vtkPolyDataWriter> pdWriter =
 vtkSmartPointer<vtkPolyDataWriter>::New();
 vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
 pd->SetPoints(pts);

 vtkSmartPointer<vtkSurfaceReconstructionFilter> surf =
 vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
 surf->SetInputData(pd);

 vtkSmartPointer<vtkContourFilter> cf =
 vtkSmartPointer<vtkContourFilter>::New();
 cf->SetInputConnection(surf->GetOutputPort());
 cf->SetValue(0, 0.0);

 // Sometimes the contouring algorithm can create a volume whose gradient
 // vector and ordering of polygon (using the right hand rule) are
 // inconsistent. vtkReverseSense cures this problem.
 vtkSmartPointer<vtkReverseSense> reverse =
 vtkSmartPointer<vtkReverseSense>::New();
 reverse->SetInputConnection(cf->GetOutputPort());
 reverse->ReverseCellsOn();
 reverse->ReverseNormalsOn();
 reverse->Update();

 pdWriter->SetInputConnection(reverse->GetOutputPort());
 pdWriter->SetFileName("PoissonSurface.vtk");
 pdWriter->Write();

 vtkSmartPointer<vtkDataSetSurfaceFilter> sf
 = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
 sf->SetInputConnection(reverse->GetOutputPort());
 sf->Update();

 vtkSmartPointer<vtkMeshQuality> mq = vtkSmartPointer<vtkMeshQuality>::New();
 mq->SetTriangleQualityMeasure(VTK_QUALITY_ASPECT_RATIO);
 mq->SetInputConnection(sf->GetOutputPort());
 mq->Update();
 double quality = mq->GetOutput()->GetFieldData()->
 GetArray("Mesh Triangle Quality")->GetComponent(0, 1);
 std::cout << "Average aspect ratio for triangles of new mesh = " << quality << "." << std::endl;
 pd = sf->GetOutput();

 std::map<int, int> nodeIndexMap;
 for (int i = 0; i < pd->GetNumberOfPoints(); i++) {
 double x[3] = { 0.0, 0.0, 0.0 };
 pd->GetPoint(i, &x[0]);
 tvmet::Vector<double, 3> temp(x[0], x[1], x[2]);
 for (int j = 0; j < nodes.size(); j++) {
 tvmet::Vector<double, 3> currNode = nodes[j]->point();
 if (tvmet::norm2(temp - currNode) < 1e-6) {
 nodeIndexMap[i] = j;
 break;
 }
 }
 }
 tvmet::Vector<int, 3> c;
 int ntri = pd->GetNumberOfCells();
 vector<tvmet::Vector<int, 3> > connectivities;
 connectivities.reserve(ntri);
 for (int i = 0; i < ntri; i++) {
 for (int a = 0; a < 3; a++) {
 c[a] = nodeIndexMap[pd->GetCell(i)->GetPointId(a)];
 }
 connectivities.push_back(c);
 }
 return connectivities;
 }
 */
/*
 This method returns mesh quality from a vector of DeformationNodes and a
 connectivity vector
 */
std::vector<double> getMeshQualityInfo(
		const std::vector<DeformationNode<3>*> &nodes,
		const std::vector<tvmet::Vector<int, 3> > connectivities) {
//Convert DeformationNodes to vtkPoints
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < nodes.size(); i++) {
		tvmet::Vector<double, 3> cp = nodes[i]->point();
		pts->InsertNextPoint(cp(0), cp(1), cp(2));
	}
//Convert connectivity information to a vtkCellArray
	vtkSmartPointer<vtkCellArray> cellArray =
			vtkSmartPointer<vtkCellArray>::New();
	vtkIdType cell[3] = { 0, 0, 0 };
	for (int i = 0; i < connectivities.size(); i++) {
		cell[0] = connectivities[i](0);
		cell[1] = connectivities[i](1);
		cell[2] = connectivities[i](2);
		cellArray->InsertNextCell(3, cell);
	}
//Make a polydata from points and cells
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->SetPoints(pts);
	pd->SetPolys(cellArray);
//Measure mesh quality -> Aspect Ratio
	vtkSmartPointer<vtkMeshQuality> mq = vtkSmartPointer<vtkMeshQuality>::New();
	mq->SetTriangleQualityMeasure(VTK_QUALITY_RADIUS_RATIO);
	mq->SetInputData(pd);
	mq->Update();
	double avgQuality = mq->GetOutput()->GetFieldData()->GetArray(
			"Mesh Triangle Quality")->GetComponent(0, 1);
	double minQuality = mq->GetOutput()->GetFieldData()->GetArray(
			"Mesh Triangle Quality")->GetComponent(0, 0);
	double maxQuality = mq->GetOutput()->GetFieldData()->GetArray(
			"Mesh Triangle Quality")->GetComponent(0, 2);
	std::vector<double> quality = { minQuality, avgQuality, maxQuality };
	return quality;
}
/*
 The following function plots MorseBonds
 */
void plotMorseBonds(const std::vector<std::string> &fileNames,
		std::string fname, double epsilon, double Rshift, double sigma,
		vtkSmartPointer<vtkCellArray> bonds) {
	std::stringstream sstm;
	std::string fileName;
	for (int z = 0; z < fileNames.size(); z++) {
		fileName = fileNames[z];
		//Check that the file exists
		assert(ifstream(fileName.c_str()));
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
				vtkPolyDataReader>::New();
		reader->SetFileName(fileName.c_str());
		vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();
		reader->Update();
		vtkSmartPointer<vtkDataArray> displacements =
				mesh->GetPointData()->GetVectors("displacements");
		vtkSmartPointer<vtkDoubleArray> force =
				vtkSmartPointer<vtkDoubleArray>::New();
		force->SetNumberOfTuples(bonds->GetNumberOfCells());
		force->SetNumberOfComponents(1);
		force->SetName("MutualForce");

		//Iterate over cells (lines) in bonds and calculate forces
		vtkSmartPointer<vtkIdList> points = vtkSmartPointer<vtkIdList>::New();
		int forceIdx = 0;
		bonds->InitTraversal();
		while (bonds->GetNextCell(points)) {
			tvmet::Vector<double, 3> r1(0.0), r2(0.0), d1(0.0), d2(0.0);
			vtkIdType a, b;
			a = points->GetId(0);
			b = points->GetId(1);
			mesh->GetPoint(a, &r1[0]);
			mesh->GetPoint(b, &r2[0]);
			displacements->GetTuple(a, &d1[0]);
			displacements->GetTuple(b, &d2[0]);
			double r = tvmet::norm2(r1 + d1 - r2 - d2);
			double forceMag = -((2.0 * sigma * exp(-(r - Rshift) * sigma)
					- 2.0 * sigma * exp(-2.0 * (r - Rshift) * sigma)) * epsilon);
			force->SetTuple1(forceIdx++, forceMag);
		}
		vtkSmartPointer<vtkPolyData> out = vtkSmartPointer<vtkPolyData>::New();
		out->SetPoints(mesh->GetPoints());
		out->SetLines(bonds);
		out->GetCellData()->AddArray(force);
		out->GetPointData()->AddArray(displacements);
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<
				vtkPolyDataWriter>::New();
		writer->SetInputData(out);
		sstm << fname << "-MorseBonds-" << z << ".vtk";
		fileName = sstm.str();
		sstm.str("");
		sstm.clear();
		writer->SetFileName(fileName.c_str());
		writer->Write();
	}
}

/*
 The next function is an overloaded version of getRadialStats(). It uses many more
 points from LoopShellSurface to calculate asphericity etc.
 */
std::vector<double> getRadialStats(vtkSmartPointer<vtkPolyData> pd,
		tvmet::Vector<double, 3> Xavg) {
	std::vector<double> qpRadius(pd->GetNumberOfPoints(), 0.0);
#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (int e = 0; e < pd->GetNumberOfPoints(); e++) {
		tvmet::Vector<double, 3> Xq(0.0);
		pd->GetPoint(e, &Xq[0]);
		double qpR = tvmet::norm2(Xq - Xavg);
		qpRadius[e] = qpR;
	}
	double Ravg = 0.0;
	for (int i = 0; i < qpRadius.size(); i++) {
		Ravg += qpRadius[i];
	}
	Ravg /= qpRadius.size();

	double dRavg2 = 0.0;
	for (int i = 0; i < qpRadius.size(); i++) {
		double dR = qpRadius[i] - Ravg;
		dRavg2 += dR * dR;
	}
	dRavg2 /= qpRadius.size();

	double asphericity = dRavg2 / (Ravg * Ravg);

	tvmet::Vector<double, 3> X0(0.0), X1(0.0);
	pd->GetPoint(0, &X0[0]);
	double closestNeighDist = 1e16;
	for (int j = 1; j < pd->GetNumberOfPoints(); j++) {
		pd->GetPoint(j, &X1[0]);
		double dist = tvmet::norm2(X1 - X0);
		if (dist < closestNeighDist) {
			closestNeighDist = dist;
		}
	}
	std::vector<double> output = { Ravg, asphericity, closestNeighDist };
	return output;
}

/*
 This function is meant to calculate theta and phi limits for bins on
 a spherical surface. It is meant to avoid common code from multiple
 drivers
 */
std::vector<vector<double> > getSphCellLimits(
		const vtkSmartPointer<vtkPolyData> &pd, int long_res) {

	bool debug = false;
	vtkSmartPointer<vtkIdList> points = vtkSmartPointer<vtkIdList>::New();
	std::vector<vector<double> > cellLimits;
	vtkSmartPointer<vtkCellArray> bins = pd->GetPolys();
	bins->InitTraversal();
	int cellId = 0;

	while (bins->GetNextCell(points)) {

		if (debug) {
			std::cout << "Cell Id: " << cellId++ << std::endl;
		}

		double theta_max = 0, theta_min = 180, phi_max = 0, phi_min = 360;

		for (int i = 0; i < points->GetNumberOfIds(); i++) {

			if (debug) {
				std::cout << "\tPoint Id : " << points->GetId(i) << std::endl
						<< "\t\t";
			}

			//Get Cartesian coordinates for each point
			double *xyz = pd->GetPoint(points->GetId(i));

			if (debug) {
				std::cout << xyz[0] << "," << xyz[1] << "," << xyz[2]
						<< std::endl << "\t\t";
			}

			//Skip the "poles" of the sphere
			if ((std::abs(xyz[0]) < 1e-8) && (std::abs(xyz[1]) < 1e-8)) {
				if (std::abs(xyz[2] - 1) < 1e-8)
					theta_min = 0;
				if (std::abs(xyz[2] + 1) < 1e-8)
					theta_max = 180;
				if (debug) {
					std::cout << std::endl;
				}
				continue;
			}

			//Convert to spherical coordinates (phi,theta)
			double phi = (180 / M_PI) * atan2(xyz[1], xyz[0]);
			double theta = (180 / M_PI) * acos(xyz[2]);

			phi = (phi < 0) ? (360 + phi) : phi;

			if (debug) {
				std::cout << "Phi = " << phi << " Theta = " << theta
						<< std::endl;
			}

			//Compare to update max and min values
			theta_max = std::max(theta, theta_max);
			theta_min = std::min(theta, theta_min);
			phi_min = std::min(phi, phi_min);
			phi_max = std::max(phi, phi_max);
		}
		//Checking for the last bin along phi direction
		if ((phi_max - phi_min) > 2 * (360 / (long_res - 1.0))) {
			phi_min = phi_max;
			phi_max = 360;
		}
		vector<double> temp;
		temp.push_back(phi_min);
		temp.push_back(phi_max);
		temp.push_back(theta_min);
		temp.push_back(theta_max);

		if (debug) {
			std::cout << "\tPhi_min_max: " << phi_min << "," << phi_max
					<< std::endl << "\t" << "Theta_min_max: " << theta_min
					<< "," << theta_max << std::endl;
		}

		cellLimits.push_back(temp);

	}
	if (debug) {
		std::cout << "Printing the bins : " << std::endl;
		std::cout << "\tBinId\tPhi_min\tPhi_max\tTheta_min\tTheta_max"
				<< std::endl;
		for (int binIter = 0; binIter < cellLimits.size(); binIter++) {
			std::cout << "\t" << binIter << "\t" << cellLimits[binIter][0]
					<< "\t" << cellLimits[binIter][1] << "\t"
					<< cellLimits[binIter][2] << "\t" << cellLimits[binIter][3]
					<< std::endl;
		}
	}
	return cellLimits;
}
/*
 This function classifies the particles belonging to certain bins
 */
void putParticlesInBins(const std::vector<std::vector<double> > &cellLimits,
		const Eigen::Matrix3Xd &newCurr, int size,
		vtkSmartPointer<vtkDoubleArray> binDensity, int viterMax) {

	bool debug = false;
	for (int i = 0; i < size; i++) {

		if (debug) {
			std::cout << "\tPoint Id = " << i << std::endl;
		}

		tvmet::Vector<double, 3> pos(0.0);
		for (int row = 0; row < 3; row++) {
			pos(row) = newCurr(row, i);
		}

		if (debug) {
			std::cout << "\t\tOriginal : " << pos << std::endl;
		}

		tvmet::Vector<double, 3> normalizedPos(0.0);
		normalizedPos = pos / tvmet::norm2(pos);

		if (debug) {
			std::cout << "\t\tNormalized : " << normalizedPos << std::endl;
		}

		//Convert to spherical coordinates (phi,theta)
		double phi = (180 / M_PI) * atan2(normalizedPos(1), normalizedPos(0));
		double theta = (180 / M_PI) * acos(normalizedPos(2));

		phi = (phi < 0) ? (360 + phi) : phi;

		if (debug) {
			std::cout << "\t\tPhi = " << phi << " Theta = " << theta
					<< std::endl;
		}

		for (int binId = 0; binId < cellLimits.size(); binId++) {
			double p_min, p_max, t_min, t_max;
			p_min = cellLimits[binId][0];
			p_max = cellLimits[binId][1];
			t_min = cellLimits[binId][2];
			t_max = cellLimits[binId][3];

			if ((p_min <= phi && phi < p_max)
					&& (t_min <= theta && theta < t_max)) {

				if (debug) {
					std::cout << "\t\tBin found : " << binId << std::endl;
				}

				double tempCount = binDensity->GetTuple1(binId);
				//Size of fileNames vector corresponds to number of time steps
				binDensity->SetTuple1(binId, tempCount + (1.0 / viterMax));
				break;
			}
		}
	}
}

/*
 This function classifies the particles belonging to certain bins
 */
void putParticlesInBins(const std::vector<std::vector<double> > &cellLimits,
		const Matrix6Xd &newCurr, int size,
		vtkSmartPointer<vtkDoubleArray> binDensity, int viterMax) {

	bool debug = false;
	for (int i = 0; i < size; i++) {

		if (debug) {
			std::cout << "\tPoint Id = " << i << std::endl;
		}

		tvmet::Vector<double, 3> pos(0.0);
		for (int row = 0; row < 3; row++) {
			pos(row) = newCurr(row, i);
		}

		if (debug) {
			std::cout << "\t\tOriginal : " << pos << std::endl;
		}

		Vector3D normalizedPos(0.0);
		normalizedPos = pos / tvmet::norm2(pos);

		if (debug) {
			std::cout << "\t\tNormalized : " << normalizedPos << std::endl;
		}

		//Convert to spherical coordinates (phi,theta)
		double phi = (180 / M_PI) * atan2(normalizedPos(1), normalizedPos(0));
		double theta = (180 / M_PI) * acos(normalizedPos(2));

		phi = (phi < 0) ? (360 + phi) : phi;

		if (debug) {
			std::cout << "\t\tPhi = " << phi << " Theta = " << theta
					<< std::endl;
		}

		for (int binId = 0; binId < cellLimits.size(); binId++) {
			double p_min, p_max, t_min, t_max;
			p_min = cellLimits[binId][0];
			p_max = cellLimits[binId][1];
			t_min = cellLimits[binId][2];
			t_max = cellLimits[binId][3];

			if ((p_min <= phi && phi < p_max)
					&& (t_min <= theta && theta < t_max)) {

				if (debug) {
					std::cout << "\t\tBin found : " << binId << std::endl;
				}

				double tempCount = binDensity->GetTuple1(binId);
				//Size of fileNames vector corresponds to number of time steps
				binDensity->SetTuple1(binId, tempCount + (1.0 / viterMax));
				break;
			}
		}
	}
}

/*
 The method insertValenceInVtk() inserts valence information in a vtk file.
 The calling method must ensure that 'fileName' exists and is a valid vtk file.
 MUST use the extension '.vtk' in 'fileName'.
 'mesh' is a pointer to a vtkPolyData
 */
void insertValenceInVtk(std::vector<std::string> &fileNames) {

#ifdef _OPENMP
# pragma omp parallel for
#endif
	for (int i = 0; i < fileNames.size(); i++) {

		//Check that the file exists
		assert(ifstream(fileNames[i].c_str()));

		//vtkDataSetReader * reader = vtkDataSetReader::New();

		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
				vtkPolyDataReader>::New();
		reader->SetFileName(fileNames[i].c_str());
		reader->ReadAllScalarsOn();
		reader->ReadAllVectorsOn();

		vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
		reader->Update();

		int numPoints = pd->GetNumberOfPoints();

		//The following vtkUnsignedIntArray will be used to store the
		//number of CELLS in the mesh that share the POINT denoted by the
		//index of the vector
		vtkSmartPointer<vtkUnsignedIntArray> countPointCells = vtkSmartPointer<
				vtkUnsignedIntArray>::New();
		countPointCells->SetNumberOfComponents(1);
		countPointCells->SetNumberOfTuples(numPoints);
		countPointCells->SetName("Valence");

		//cellIds will be used to temporarily hold the CELLS that use a
		//point specified by a point id.
		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

		pd->BuildLinks();
		for (int p = 0; p < pd->GetNumberOfPoints(); p++) {
			pd->GetPointCells(p, cellIds);
			countPointCells->SetValue(p, cellIds->GetNumberOfIds());
			cellIds->Reset();
		}

		pd->GetPointData()->AddArray(countPointCells);
		//We will use a vtkPolyDataWriter to write our modified output
		//files that will have Capsomer information as well
		vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<
				vtkPolyDataWriter>::New();
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

/*
 The following function creates a Morse bond network based on a
 search radius from a cloud of DeformationNodes. The output is
 returned in the input vtkCellArray pointer.
 */
void getMorseBonds(vtkSmartPointer<vtkCellArray> bonds,
		std::vector<DeformationNode<3>*> morseNodes, double searchRad) {
	vtkIdType bond[2] = { 0, 0 };
	for (int i = 0; i < morseNodes.size(); i++) {
		tvmet::Vector<double, 3> centerNode = morseNodes[i]->point();
		for (int a = i + 1; a < morseNodes.size(); a++) {
			tvmet::Vector<double, 3> currNode = morseNodes[a]->point();
			double r = tvmet::norm2(centerNode - currNode);
			if (r <= searchRad) {
				bond[0] = i;
				bond[1] = a;
				bonds->InsertNextCell(2, bond);
			}
		}
	}
}

}
