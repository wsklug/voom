// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Andrew R. Missel
//                University of California Los Angeles
//                 (C) 2009 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__IntersectionFinder_h__)
#define __IntersectionFinder_h__

#include "Node.h"
#include "Element.h"
#include "VoomMath.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
namespace voom
{

  template<int N>
  class IntersectionFinder : public Element {
    
  public: 

    typedef tvmet::Vector<double,N> VectorND;
    typedef DeformationNode<N> Node_t;    

    IntersectionFinder() {}

    ~IntersectionFinder() {}
    
    // determine whether two line segments intersect, return true or false and store intersection point //
    static bool checkIntersection(VectorND & p1A, VectorND & p1B, VectorND & p2A, VectorND & p2B, VectorND & intersectPt) {
      double a1,a2;
      if(abs(p1B[0]-p1A[0]) > 1.0e-6 && abs(p2B[0]-p2A[0]) > 1.0e-6) {
	a1 = (p1B[1]-p1A[1])/(p1B[0]-p1A[0]);
	a2 = (p2B[1]-p2A[1])/(p2B[0]-p2A[0]);
	intersectPt[0] = (p1A[1]-p2A[1]+a2*p2A[0]-a1*p1A[0])/(a2-a1);
	intersectPt[1] = (a2*(p1A[1]-a1*p1A[0])-a1*(p2A[1]-a2*p2A[0]))/(a2-a1);
      }
      else if(abs(p1B[0]-p1A[0]) <= 1.0e-6 && abs(p2B[0]-p2A[0]) > 1.0e-6){
	intersectPt[0] = p1A[0];
	a2 = (p2B[1]-p2A[1])/(p2B[0]-p2A[0]);
	intersectPt[1] = a2*p1A[0] + p2A[1] - a2*p2A[0];
      }
      else if(abs(p2B[0]-p2A[0]) <= 1.0e-6 && abs(p1B[0]-p1A[0]) > 1.0e-6) {
	intersectPt[0] = p2A[0];
	a1 = (p1B[1]-p1A[1])/(p1B[0]-p1A[0]);
	intersectPt[1] = p1A[1] - a1*(p1A[0] - p2A[0]);
      }
      else {
	intersectPt[0] = (p1A[0]+p2A[0])/2.0;
	intersectPt[1] = 1.0e30;
      }

      if((intersectPt[0] >= min(p1A[0],p1B[0])) && (intersectPt[0] <=  max(p1A[0],p1B[0])) && (intersectPt[1] >= min(p1A[1],p1B[1])) && (intersectPt[1] <= max(p1A[1],p1B[1])) && (intersectPt[0] >= min(p2A[0],p2B[0])) && (intersectPt[0] <=  max(p2A[0],p2B[0])) && (intersectPt[1] >= min(p2A[1],p2B[1])) && (intersectPt[1] <= max(p2A[1],p2B[1]))) {

	return true;
      }
      
      else {
	return false;
      }
    }
    
    // determine whether two line segments are close to one another, return true or false //
    static bool checkIntersection(VectorND & p1A, VectorND & p1B, VectorND & p2A, VectorND & p2B, VectorND & intersectPt1, VectorND & intersectPt2, double tol, PeriodicBox * bx) {
      // check if the segments intersect; return true if they are //
      // deprecated for now!! //

      // if segments do not intersect, check whether midpoints are within tol of each other; return true if they are //
      
      VectorND midpoint1;
      VectorND midpoint2;
      midpoint1 = (p1A+p1B)/2.0;
      midpoint2 = (p2A+p2B)/2.0;
      VectorND sep;
      sep = midpoint2 - midpoint1;
      bx->mapDistance(sep);
      if(norm2(sep) < tol) {
	intersectPt1 = midpoint1;
	intersectPt2 = midpoint2;
	return true;
      }
      else {
	return false;
      }
    }

  };
};

#endif // __IntersectionFinder_h__
