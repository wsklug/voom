// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.5  2005/05/23 17:43:20  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file ShellGeometry.h

  \brief ShellGeometry is a structure which computes and contains all
  of the geometric quantities relevant to calculations on a curved surface.

*/

#if !defined(__ShellGeometry_h__)
#define __ShellGeometry_h__

#include<blitz/array.h>
 
 
#include<vector>
#include<iostream>
#include "voom.h"

namespace voom
{

  
  /*!  Structure which computes and contains all of the geometric
    quantities relevant to calculations on a curved surface.
  */
  class ShellGeometry 
  {
  private:
    tvmet::Vector< Vector3D, 2 >    _a;
    tvmet::Matrix< Vector3D, 2, 2 > _aPartials;
    Vector3D   	                    _d;
    tvmet::Vector< Vector3D, 2 >    _dPartials;
    tvmet::Vector< Vector3D, 2 >    _aDual;
    tvmet::Matrix< Vector3D, 2, 2 > _aDualPartials;
    
    Tensor2D _metricTensor;
    Tensor2D _metricTensorInverse;

    double _metric;
    double _metricInverse;
		
  public:
    // Default Constructor
    ShellGeometry() {
      Vector3D zero;
      zero = 0.0, 0.0, 0.0;
      _a = tvmet::Vector< Vector3D, 2>(zero, zero);
      _aPartials = zero, zero, zero, zero;
      _d = zero;
      _dPartials = tvmet::Vector< Vector3D, 2>(zero, zero);
      _aDual = tvmet::Vector< Vector3D, 2>(zero, zero);
      _aDualPartials = zero, zero,zero, zero;
      _metricTensor = 
	0.0, 0.0,
	0.0, 0.0;
      _metricTensorInverse = 
	0.0, 0.0,
	0.0, 0.0;
      _metric = 0.0;
      _metricInverse = 0.0;
    };
    
    //! Constructor
    ShellGeometry( const tvmet::Vector< Vector3D, 2 > &	a, 
		   const tvmet::Matrix< Vector3D, 2, 2 > & aPartials ) {
      update( a, aPartials );
    }
    
    void update( const tvmet::Vector< Vector3D, 2 >& 	a,
		 const tvmet::Matrix< Vector3D, 2, 2 >& aPartials ) {
      // recompute everything from a and aPartials 
    
      _a = a;
      _aPartials = aPartials;

      _updateMetrics();
      _updateDual();
      _updateDirector();
      _updateDirectorPartials();
      _updateDualPartials();			
    }

    const tvmet::Vector<Vector3D,2>&   a()	   const {return _a;}
    const tvmet::Matrix<Vector3D,2,2>& aPartials() const {return _aPartials;}
    const Vector3D&		       d()	   const {return _d;}
    const tvmet::Vector<Vector3D,2>&   dPartials() const {return _dPartials;}
    const tvmet::Vector<Vector3D,2>&   aDual()	   const {return _aDual;}
    const tvmet::Matrix<Vector3D,2,2>& aDualPartials() const {return _aDualPartials;}
    
    const Tensor2D& metricTensor()	  const {return _metricTensor;}
    const Tensor2D& metricTensorInverse() const {return _metricTensorInverse;}

    double metric()	const  {return _metric;}
    double metricInverse() const {return _metricInverse;}
		
    //! overload operator =
    ShellGeometry& operator = ( const ShellGeometry& g );
	  
  private:
    void _updateMetrics();
    void _updateDual();
    void _updateDirector();
    void _updateDualPartials();
    void _updateDirectorPartials();

  };


} // namespace voom
#endif // __ShellGeometry_h__
