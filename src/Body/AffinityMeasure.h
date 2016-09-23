// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Mo Bai
//                University of California Los Angeles
//                   (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__AffinityMeasure_h__)
#define __AffinityMeasure_h__

#include <blitz/array.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "voom.h"
#include "Node.h"

#include "vtkPolyData.h"
#include "vtkDelaunay2D.h"
#include "vtkCellArray.h"
#include "vtkShrinkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRegressionTestImage.h"
#include <vtkPolyDataWriter.h>
#include "AffinityElement.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

namespace voom
{

  class AffinityMeasure
  {
  public:
    
    typedef std::vector<Tensor2D> StrainContainer;
    typedef std::vector<Vector2D> CentroidContainer;
    typedef BrownianNode<2> Node_t;
    typedef std::vector<Node_t *> NodeContainer;
    typedef NodeContainer::iterator NodeIterator;
    typedef NodeContainer::const_iterator ConstNodeIterator;
    typedef std::vector<AffinityElement *> AffinityElementContainer;
    typedef std::vector< std::pair<Vector2D,Tensor2D> > StrainField;


    //! Construct from stuff
    AffinityMeasure(const NodeContainer & nodes);
    
    //! virtual destructor
    virtual ~AffinityMeasure() {;}

    int triangulate();

    //get number of elements
    int getnElements () const;

    //StrainContainer getStrains() const {return _strains;}

    //CentroidContainer getCentroids() const {return _centroids;}

    AffinityElementContainer & getElements();

    void resetNodes(NodeContainer & nodes);

    StrainField getStrainField();

    double strainMeasure(double affShear, double maxArea);

    double rotationMeasure(double affRotation, double maxArea);

    void printParaview(const std::string name) const;

  private:
    NodeContainer _nodes;
    StrainContainer _strains;
    CentroidContainer _centroids;
    int _nElements;
    AffinityElementContainer _elements;
    std::vector< tvmet::Vector<int,3> > _connectivities;
  };  

} // namespace voom

//#include "AffinityMeasure.icc"

#endif // __AffinityMeasure_h__
