// -*- C++ -*-
#include <vector>
#include <tvmet/Vector.h>
#include "TriangleQuadrature.h"
#include "LinearizedElement2D.h"
#include "Hookean2D.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "Node.h"

using namespace std;
using namespace voom;

int main()
{
  DeformationNode<2>::PositionVector x3[3];
  DeformationNode<2>::Point u3[3];

  x3[0] = -0.50, 0.00;	  u3[0] = -0.25, 0.00;
  x3[1] =  0.50, 0.00;	  u3[1] =  0.25, 0.00;
  x3[2] =  0.00, 1.00;	  u3[2] =  0.00, 0.00;

  typedef LinearizedElement2D<  
    DeformationNode<2>, 
    TriangleQuadrature, 
    Hookean2D, 
    ShapeTri3          > Tri3Element;

  Tri3Element::NodeContainer nodes3;
  for(int i=0; i<3; i++) {
    NodeBase::DofIndexMap idx(2);
    idx[0] = 2*i; idx[1] = 2*i+1;
    DeformationNode<2> * n = new DeformationNode<2>( i, idx, x3[i], u3[i] );
    nodes3.push_back( n );
  }

  TriangleQuadrature quad3(1);
  Hookean2D mat(Hookean2D::planeStrain, 10.0, 0.3);

  // global index map
  Tri3Element::DofIndexMap index3(6);
  for(int j=0; j<6; j++) index3[j] = j;   
  Tri3Element tri3(quad3,mat,nodes3);

  tri3.compute(true,true,true);

  for( Tri3Element::ConstQuadPointIterator p=tri3.quadraturePoints().begin(); 
       p!=tri3.quadraturePoints().end(); p++ ) {
    std::cout <<"Strain = "<< p->material.strain() << std::endl
	      <<"Stress = "<< p->material.stress() << std::endl;
  }

  tri3.checkConsistency();
  tri3.checkRank(3);


  DeformationNode<2>::PositionVector x6[6];
  DeformationNode<2>::Point u6[6];

//   x6[0] =  0.00, 0.00;	  u6[0] =  0.00, 0.00;
//   x6[1] =  1.00, 0.00;	  u6[1] =  0.50, 0.00;
//   x6[2] =  0.00, 1.00;	  u6[2] =  0.00, 0.00;
//   x6[3] =  0.50, 0.00;	  u6[3] =  0.25, 0.00;
//   x6[4] =  0.50, 0.50;	  u6[4] =  0.25, 0.00;
//   x6[5] =  0.00, 0.50;	  u6[5] =  0.00, 0.00;
  x6[0] =  0.00, 0.00;	  u6[0] =  0.0000, 0.00;
  x6[1] =  0.50, 0.50;	  u6[1] =  0.0050, 0.00;
  x6[2] =  0.00, 1.00;	  u6[2] =  0.0000, 0.00;
  x6[3] =  0.25, 0.25;	  u6[3] =  0.0025, 0.00;
  x6[4] =  0.25, 0.75;	  u6[4] =  0.0025, 0.00;
  x6[5] =  0.00, 0.50;	  u6[5] =  0.0000, 0.00;

  typedef LinearizedElement2D<  
    DeformationNode<2>, 
    TriangleQuadrature, 
    Hookean2D, 
    ShapeTri6          > Tri6Element;

  Tri6Element::NodeContainer nodes6;
  for(int i=0; i<6; i++) {
    NodeBase::DofIndexMap idx(2);
    idx[0] = 2*i; idx[1] = 2*i+1;
    DeformationNode<2> * n = new DeformationNode<2>( i, idx, x6[i], u6[i] );
    nodes6.push_back( n );
  }

  TriangleQuadrature quad6(2);

  // global index map
  Tri3Element::DofIndexMap index6(12);
  for(int j=0; j<12; j++) index6[j] = j;   
  Tri6Element tri6(quad6,mat,nodes6);

  tri6.compute(true,true,true);

  for( Tri6Element::ConstQuadPointIterator p=tri6.quadraturePoints().begin(); 
       p!=tri6.quadraturePoints().end(); p++ ) {
    std::cout <<"Strain = "<< p->material.strain() << std::endl
	      <<"Stress = "<< p->material.stress() << std::endl;
  }

  tri6.checkConsistency();
  tri6.checkRank(9);

  for(vector< DeformationNode<2>* >::iterator n=nodes3.begin(); 
      n != nodes3.end(); n++ ) delete *n;

  for(vector< DeformationNode<2>* >::iterator n=nodes6.begin(); 
      n != nodes6.end(); n++ ) delete *n;

  return 0;
}
