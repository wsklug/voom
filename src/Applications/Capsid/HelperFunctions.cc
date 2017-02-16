#include "HelperFunctions.h"

/*
	The method writeEdgeStrainVtk() inserts edge strain information in a vtk file
	The calling method must ensure that 'fileName' exists and is a valid vtk file.
	MUST use the extension '.vtk' in 'fileName'.
*/
namespace voom
{
	using namespace std;

	void writeEdgeStrainVtk(std::vector<std::string> &fileNames,
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

	void writeEdgeStrainVtk(std::vector<std::string> &fileNames,
		double avgEdgeLen, std::vector<double> percentStrain) {

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

	std::vector<double> calcEdgeLenAndStdDev
	(const std::vector< DeformationNode<3>* > &defNodes,
		const vector< tvmet::Vector<int, 3> > &connectivities) {

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
	/*
		This function returns a surface made of a cloud of points
	*/
	vector<tvmet::Vector<int, 3> > delaunay3DSurf(const std::vector<DeformationNode<3>*> 
		&nodes)
	{
		vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
		for (int i = 0; i < nodes.size(); i++) {
			tvmet::Vector<double, 3> cp = nodes[i]->point();
			pts->InsertNextPoint(cp(0), cp(1), cp(2));
		}
		vtkSmartPointer<vtkPolyDataWriter> pdWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
		vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
		pd->SetPoints(pts);

		vtkSmartPointer<vtkDelaunay3D> d3D = vtkSmartPointer<vtkDelaunay3D>::New();
		d3D->SetInputData(pd);
		vtkSmartPointer<vtkDataSetSurfaceFilter> sf = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
		sf->SetInputConnection(d3D->GetOutputPort());
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

	/*
		This method creates a 3D surface from unorganized points using 
		Poisson surface reconstruction. It returns a connectivity matrix
	*/
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

	/*
		This function is meant to calculate theta and phi limits for bins on
		a spherical surface. It is meant to avoid common code from multiple
		drivers
	*/
	std::vector<vector<double> > getSphCellLimits(const vtkSmartPointer<vtkPolyData> &pd, int long_res) {

		bool debug = false;
		vtkSmartPointer<vtkIdList> points
			= vtkSmartPointer<vtkIdList>::New();
		std::vector<vector<double> > cellLimits;
		vtkSmartPointer<vtkCellArray> bins = pd->GetPolys();
		bins->InitTraversal();
		int cellId = 0;

		while (bins->GetNextCell(points)) {

			if (debug) {
				std::cout << "Cell Id: " << cellId++
					<< std::endl;
			}

			double theta_max = 0, theta_min = 180,
				phi_max = 0, phi_min = 360;

			for (int i = 0; i < points->GetNumberOfIds(); i++) {

				if (debug) {
					std::cout << "\tPoint Id : " << points->GetId(i)
						<< std::endl << "\t\t";
				}

				//Get Cartesian coordinates for each point
				double *xyz = pd->GetPoint(points->GetId(i));

				if (debug) {
					std::cout << xyz[0] << "," << xyz[1] << ","
						<< xyz[2] << std::endl << "\t\t";
				}

				//Skip the "poles" of the sphere
				if ((std::abs(xyz[0]) < 1e-8) &&
					(std::abs(xyz[1]) < 1e-8)) {
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
				double phi = (180 / M_PI)*atan2(xyz[1], xyz[0]);
				double theta = (180 / M_PI)*acos(xyz[2]);

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
					<< std::endl << "\t"
					<< "Theta_min_max: " << theta_min << "," << theta_max
					<< std::endl;
			}

			cellLimits.push_back(temp);

		}
		if (debug) {
			std::cout << "Printing the bins : " << std::endl;
			std::cout << "\tBinId\tPhi_min\tPhi_max\tTheta_min\tTheta_max" << std::endl;
			for (int binIter = 0; binIter < cellLimits.size(); binIter++) {
				std::cout << "\t" << binIter
					<< "\t" << cellLimits[binIter][0]
					<< "\t" << cellLimits[binIter][1]
					<< "\t" << cellLimits[binIter][2]
					<< "\t" << cellLimits[binIter][3]
					<< std::endl;
			}
		}
		return cellLimits;
	}
	/*
		This function classifies the particles belonging to certain bins
	*/
	void putParticlesInBins(const std::vector<std::vector<double> > &cellLimits,
		const Eigen::Matrix3Xd &newCurr, const std::vector<DeformationNode<3>*> &defNodes,
		vtkSmartPointer<vtkDoubleArray> binDensity, int viterMax) {

		bool debug = false;
		for (int i = 0; i < defNodes.size(); i++) {

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
			double phi = (180 / M_PI)*atan2(normalizedPos(1),
				normalizedPos(0));
			double theta = (180 / M_PI)*acos(normalizedPos(2));

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

				if ((p_min <= phi && phi < p_max) &&
					(t_min <= theta && theta < t_max))
				{

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
