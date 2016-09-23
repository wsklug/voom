// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         Ankush Aggarwal
//                University of California Los Angeles
//                   (C) 2011 All Rights Reserved
//
//----------------------------------------------------------------------
// 

#if !defined(__RadialForce_h__)
#define __RadialForce_h__

#include <vector>
#include <cstdio>
#include <ctime>

#include "Element.h"

namespace voom {
  //! Applies a potential that provides a radial expansion/compression for spherical shaped objects
  // E = k/2 (Rmean -R0)^2

  class RadialForce :
    public Element {
   
    public:
    using Element::BaseNodeContainer; 
    using Element::BaseNodeIterator; 
    using Element::ConstBaseNodeIterator; 

    //! Construct from a set of nodes and spring coefficient k
    RadialForce( const BaseNodeContainer & nodes, double k)
      {
        _baseNodes = nodes;
        _k = k;
        center();
        _radius = MeanRadius();
      }

    void setR(double R){_radius=R;}

    void center(){
      BaseNodeIterator n=_baseNodes.begin();
      _center.resize((*n)->dof(),0.);
      for(; n!=_baseNodes.end(); n++){
	for(int i=0; i<(*n)->dof(); i++) 
          _center[i] += (*n)->getPoint(i);
      }
      for(int i=0; i<_center.size(); i++) 
	_center[i] /= _baseNodes.size();
    }
      

    double MeanRadius() {
      double meanR=0;
      center();
      for(BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
        double normx=0;
	for(int i=0; i<(*n)->dof(); i++) 
          normx += ((*n)->getPoint(i)-_center[i])*((*n)->getPoint(i)-_center[i]);
        normx = sqrt(normx);
        meanR += normx;
      }
      meanR /= _baseNodes.size();
      return meanR;
      }

    //! Compute energy, force and stiffness
    void compute(bool f0, bool f1, bool f2) {
      center();
      double Rmean=MeanRadius();
      double Rdiff = Rmean-_radius;
      
      if(f0) _energy = _k/2.*Rdiff*Rdiff;

      if(f1){
        double force=_k*Rdiff/_baseNodes.size(); //this magnitude of force is applied at each node in radial direction
        BaseNodeIterator n=_baseNodes.begin();
        std::vector<double> comm_force((*n)->dof(),0.);  //comm_force is the force added to all nodes to make sure that the total force on all nodes is zero
        for(; n!=_baseNodes.end(); n++){
          double normx=0;
	  for(int i=0; i<(*n)->dof(); i++) 
            normx += ((*n)->getPoint(i)-_center[i])*((*n)->getPoint(i)-_center[i]);
          normx = sqrt(normx);
          if(normx>1e-5){
	    for(int i=0; i<(*n)->dof(); i++) {
	      double xi = (*n)->getPoint(i)-_center[i];
	      (*n)->addForce(i, force*xi/normx);
              comm_force[i] += force*xi/normx;
	    }
          }
        }
        for(int i=0; i<comm_force.size(); i++) comm_force[i] /= _baseNodes.size(); 
        for(n=_baseNodes.begin(); n!=_baseNodes.end(); n++){
          for(int i=0; i<(*n)->dof(); i++) 
            (*n)->addForce(i, -comm_force[i]);
        }
      }
    }
   
   private:

   //! spring constant
   double _k;

   //! radius R0 (energy=0 when MeanRadius=_radius)
   double _radius;

   //! center of the nodes set to apply radial force (in case the given set is not symmetric)? - right now assumed at origin
   std::vector<double> _center;
   };

} // end namespace

#endif // __RadialForce_h__
