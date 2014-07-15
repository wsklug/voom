// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Mo Bai
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

//#include "AffinityMeasure.h"

namespace voom
{

  AffinityMeasure::AffinityMeasure(const NodeContainer & nodes){
    _nodes = nodes;
    _strains.clear();
    _centroids.clear();
    _nElements = 0;
    _elements.clear();
    _connectivities.clear();
  }

  int AffinityMeasure::triangulate(){
    vtkPoints *newPts = vtkPoints::New();
   
    for (int i=0;i<_nodes.size();i++){
      const Vector2D & x = _nodes[i]->point();
      newPts->InsertNextPoint(  x(0),  x(1), 0.0 );
    }

    vtkIdType inNumPts = newPts->GetNumberOfPoints();
    cout << "input numPts= " << inNumPts << endl;

    vtkPolyData *pointCloud = vtkPolyData::New();
    pointCloud->SetPoints(newPts);
    newPts->Delete();

    vtkDelaunay2D *delaunay2D = vtkDelaunay2D::New();
    delaunay2D->SetInput( pointCloud );
    pointCloud->Delete();
    delaunay2D->Update();

    vtkPolyData *triangulation = delaunay2D->GetOutput();

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput( triangulation );
    writer->SetFileName( "mesh.vtk" ); 
    writer->SetFileTypeToASCII();
    writer->Write();

    vtkIdType outNumPts = triangulation->GetNumberOfPoints();
    vtkIdType outNumCells = triangulation->GetNumberOfCells();
    vtkIdType outNumPolys = triangulation->GetNumberOfPolys();
    vtkIdType outNumLines = triangulation->GetNumberOfLines();
    vtkIdType outNumVerts = triangulation->GetNumberOfVerts();

    cout << "output numPts= " << outNumPts << endl;
    cout << "output numCells= " << outNumCells << endl;
    cout << "output numPolys= " << outNumPolys << endl;
    cout << "output numLines= " << outNumLines << endl;
    cout << "output numVerts= " << outNumVerts << endl;

    if( outNumPts != inNumPts )
      {
	cout << "ERROR: output numPts " << outNumPts
	     << " doesn't match input numPts=" << inNumPts << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    if( !outNumCells )
      {
	cout << "ERROR: output numCells= " << outNumCells << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    if( outNumPolys != outNumCells )
      {
	cout << "ERROR: output numPolys= " << outNumPolys
	     << " doesn't match output numCells= " << outNumCells << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    if( outNumLines )
      {
	cout << "ERROR: output numLines= " << outNumLines << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    if( outNumVerts )
      {
	cout << "ERROR: output numVerts= " << outNumVerts << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    // check that every point is connected
    triangulation->BuildLinks();

    vtkIdList *cellIds = vtkIdList::New();
    vtkIdType numUnconnectedPts = 0;

    for(vtkIdType ptId=0; ptId<outNumPts; ptId++)
      {
	triangulation->GetPointCells(ptId,cellIds);
	if( !cellIds->GetNumberOfIds() )
	  {
	    numUnconnectedPts++;
	  }
      }

    cellIds->Delete();

    cout << "Triangulation has " << numUnconnectedPts
	 << " unconnected points" << endl;

    if( numUnconnectedPts )
      {
	cout << "ERROR: Triangulation has " << numUnconnectedPts
	     << " unconnected points" << endl;
	delaunay2D->Delete();
	return EXIT_FAILURE;
      }

    tvmet::Vector<double, 2> coords;
    coords = 0.3 , 0.3; //can be any values, because they don't change the strain tensor
    ShapeTri3 s(coords);
    _nElements = triangulation->GetNumberOfCells();
    int nodes_per_element = 3; //triangle
    _strains.clear();
    _centroids.clear();

    _connectivities.clear();
    for(int e=0; e<_nElements; e++) {
      NodeContainer elementnodes;
      elementnodes.clear();
      Vector2D centroid(0.0);

      tvmet::Vector<int,3> myNodes;
      for(int a=0; a<nodes_per_element; a++) {
	int A = triangulation->GetCell(e)->GetPointId(a);
	myNodes(a)=A;
	//elementnodes.push_back(_nodes[A]);
	centroid += _nodes[A]->point();
      }
      
      // check to make sure triangle nodes are in counter-clockwise order
      Vector2D node0 = _nodes[myNodes(0)]->point();
      Vector2D node1 = _nodes[myNodes(1)]->point();
      Vector2D node2 = _nodes[myNodes(2)]->point();

      tvmet::Vector<double,2> AB, AC;
      AB = node1 - node0;
      AC = node2 - node0;
      double area = (AB(0)*AC(1) - AB(1)*AC(0))*0.5;
      if (area < 0.0){
	int tempID = myNodes(0);
	myNodes(0) = myNodes(1);
	myNodes(1) = tempID;
      }

      for(int a=0; a<nodes_per_element; a++) {
	int nodeID = myNodes(a);
	elementnodes.push_back(_nodes[nodeID]);
      }

      _connectivities.push_back(myNodes);

      centroid = centroid/nodes_per_element;
      _centroids.push_back(centroid);

      AffinityElement * element = new AffinityElement(s,elementnodes);
      _elements.push_back(element);
      Tensor2D strain(0.0);
      element->compute();
      strain = element->Strain();
      _strains.push_back(strain);
    }
    
    return 0;
  }


  AffinityMeasure::StrainField AffinityMeasure::getStrainField() {
    int nStrains = _strains.size();
    int nCentroids = _centroids.size();
    assert(nStrains == nCentroids);
    StrainField sf(nStrains);
    for(int i=0; i<nStrains; i++) {
      sf[i] = std::pair<Vector2D,Tensor2D>(_centroids[i],_strains[i]);
    }
    return sf;
  }

  double AffinityMeasure::Measure(double uniShear) const {
    double area = 0.0;
    double dShear = 0.0;
    double nonaffinity = 0.0;
    for (int i=0; i<_elements.size(); i++) {
      double dShear = _elements[i]->Rotation() - uniShear;
      double elementArea = _elements[i]->Area();
      nonaffinity += dShear*dShear*elementArea;
      area += elementArea;
    }
    nonaffinity /= area;
    return nonaffinity;
  }

  //! create input file used by Paraview, a 3D viewer
  void AffinityMeasure::printParaview(const std::string name) const
  {
    std::string fileName = name + ".vtk";
    std::ofstream ofs(fileName.c_str());

    if (!ofs) {
      std::cout << "can not open output ("
		<< fileName
		<< ") file." << std::endl;
      exit(0);
    }
    
    ////////////////////////////////////////////////////////////////////
    //
    //    Node Section
    //
    ofs << "# vtk DataFile Version 2.0\n"
	<< "Test example" << std::endl
	<< "ASCII" << std::endl
	<< "DATASET POLYDATA" << std::endl
	<< "POINTS  " << _nodes.size() << "  double" << std::endl;
    
    //
    // output nodal postions
    for(int a=0; a<_nodes.size(); a++) {
      const Vector2D & nodalPos =  _nodes[a]->position();
      ofs << std::setprecision(16) 
	  << nodalPos(0) << "  "
	  << nodalPos(1) << "  "
	  << 0.0 << std::endl;
    }

    /////////////////////////////////////////////////////////////////////
    //
    //    Element Section
    //
    ofs << "POLYGONS  " << _connectivities.size() << "  "
	<< 4*_connectivities.size() << std::endl;
    for(int e=0; e<_connectivities.size(); e++) {
      ofs << 3 << "  "
	  << std::setw(10) << _connectivities[e](0)
	  << std::setw(10) << _connectivities[e](1)
	  << std::setw(10) << _connectivities[e](2)
	  << std::endl;
    }


     	
    //////////////////////////////////////////////////////////////////////////
    //
    //  output color for each element ( corresponding to mean
    //  curvature, energy...)
    //
    ofs << "CELL_DATA    " << _connectivities.size() << std::endl;
    //
    // output for shear strain
    ofs << "SCALARS    shear_strain    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<_strains.size(); e++) {
      ofs << _strains[e](0,1) << std::endl;
    }
    ofs << std::endl;

    // output for normal xx strain
    ofs << "SCALARS    xx_strain    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<_strains.size(); e++) {
      ofs << _strains[e](0,0) << std::endl;
    }
    ofs << std::endl;

    // output for normal yy strain
    ofs << "SCALARS    yy_strain    double    1" << std::endl;
    ofs << "LOOKUP_TABLE default" << std::endl;
    for ( int e = 0; e<_strains.size(); e++) {
      ofs << _strains[e](1,1) << std::endl;
    }
    ofs << std::endl;

    ofs << endl << "POINT_DATA " << _nodes.size() << endl
	<< "VECTORS displacements double" << endl;
    /*     << "LOOKUP_TABLE displacements" << endl; */
    //
    // output nodal postions
    for(int a=0; a<_nodes.size(); a++) {
      Vector2D nodalDisp;
      nodalDisp = _nodes[a]->point() - _nodes[a]->position();
      ofs << std::setprecision(16) 
	  << nodalDisp(0)
	  << '\t' <<nodalDisp(1)
	  << '\t' <<0.0 << std::endl;
    }
    
    ofs.close();

    return;
  }

} // namespace voom
