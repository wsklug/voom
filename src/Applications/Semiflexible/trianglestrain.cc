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

#include "Node.h"
#include "SemiflexibleGel.h"
#include <random/uniform.h>
//#include "AffinityElement.h"
#include "AffinityMeasure.h"

int main( int argc, char* argv[] ) {
  AffinityMeasure::NodeContainer nodes;
  nodes.clear();

  int a;
  unsigned id;
  BrownianNode<2>::Point X(0.0),x(0.0);
  NodeBase::DofIndexMap idx(2);
  SemiflexibleGel<2>::DefNode * nd;
  double L = 10.0;

  ranlib::Uniform<double> rnu;
  int nNodes = 12;
  for ( int a=0; a<nNodes; a++) {
    rnu.seed((unsigned int)time(0)+a);
    id = a;
    idx[0]=2*a; idx[1]=2*a+1;
    X = rnu.random()*L, rnu.random()*L;
    x = X + X * 0.01 * rnu.random();
    nd = new SemiflexibleGel<2>::DefNode(id,idx,X,x);
    nd -> setId(id);
    nodes.push_back( nd );
    std::cout << "node " << a << "= " << nd->point() << std::endl;
  }
  
  AffinityMeasure * measure = new AffinityMeasure(nodes);
  measure->triangulate();
  
  int nElements = measure->getnElements();
  std::cout << "Number of elements: " << nElements << std::endl;
  std::vector<Tensor2D> strains;
  strains = measure->getStrains();
  std::cout << "Number of strains: " << strains.size() << std::endl;
  std::vector<Vector2D> centroids;
  centroids = measure->getCentroids();
  std::cout << "Number of centroids: " << centroids.size() << std::endl;
  
  for (int i=0; i< nElements; i++){
    std::cout << "Element " << i << ":" << std::endl;
    std::cout << "Strain: " << strains[i] << std::endl;
    std::cout << "Centroid: " << centroids[i] <<  std::endl;
  }  

  return 0;
}
