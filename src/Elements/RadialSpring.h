// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2016 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Created 2016/08/04 19:39:00  amit
// 
//
//
//----------------------------------------------------------------------

#if !defined(__RadialSpring_h__)
#define __RadialSpring_h__

#include <blitz/array.h>
#include <vector>
#include "Node.h"
#include "Element.h"

namespace voom
{

  //! RadialSpring element adds a radial spring to a shell centered at
  //! origin.
  /*! RadialSpring element tries to keep all nodes on surface of a
      sphere of specified radius
   */
 
  class RadialSpring  : public Element
  {
    
  public:

    typedef DeformationNode<3> Node; // nickname for mechanics nodes
    typedef std::vector< Node* > NodeContainer;    

    //! Constructor
    RadialSpring( const NodeContainer &nodes, double springConst, double Radius);

    //! Destructor
    ~RadialSpring();
    
    //! Do mechanics on element; compute energy, forces, and/or stiffness.
    void compute(bool f0, bool f1, bool f2);

    //! Get the spring constant
    double getSpringConstant(){ return _k;}

    //! Get the radius of shell
    double getRadius(){ return _R; }

    //! Update the spring constant
    void updateSpringConstant(double k){ _k = k;}

    //! Update the radius of the shell
    void updateRadius(double R){ _R = R;}

    //! Access the container of nodes
    const NodeContainer& nodes() const { return _nodes; }    
    
  private:
    NodeContainer _nodes;
    double _R; //Radius of the shell
    double _k; //Spring constant
    int _nodeCount;
  };

} // namespace voom

#endif // __BrownianKick_h__
