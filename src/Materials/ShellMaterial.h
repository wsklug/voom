// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------

/*! 
  \file SCElastic.h

  \brief Base Class for a material object.  This template can be
  combined with generic 3-D materials for use with thin shells.

*/

#if !defined(__ShellMaterial_h__)
#define __ShellMaterial_h__

#include<vector>
#include "voom.h"
#include "ShellGeometry.h"

namespace voom
{

  /*!  Base Class template for a thin-shell material object.  This
    template can be combined with generic 3-D materials for use with
    thin shells.
  */

  template< class BaseMaterial_t >
  class ShellMaterial
  {
  public:
    // ShellMaterial() {
    //   _W = 0.0;
    //   _n(0) = _n(1) = _n(2) =
    // 	_nTC(0) = _nTC(1) = _nTC(2) = Vector3D(0.0);
    //   _m(1) = _m(2) = Vector3D(0.0);
    //   _stretch = 0.0;
    // }

    //! virtual destructor
    virtual ~ShellMaterial() { ; }
    //! initialize the geometry in reference coords
    void setRefGeometry( const ShellGeometry& g ) {_referenceGeometry = g;}
    void setGeometry( const ShellGeometry& g ) {_deformedGeometry = g;}


    virtual double energyDensity() const { return _W; }    

    virtual double stretchingEnergy() const = 0;   

    virtual double bendingEnergy() const = 0;

    virtual const tvmet::Vector< Vector3D, 3 >& stressResultants() const 
    { return _n; }
    virtual const tvmet::Vector< Vector3D, 2 >& momentResultants() const 
    { return _m; }

    virtual const tvmet::Vector< Vector3D, 3 >& totalCurvatureStressResultants() const 
    { return _nTC; }

    const double& stretch() const { return _stretch; }

    const ShellGeometry& shellGeometry() const { return _deformedGeometry; }
    const ShellGeometry& refShellGeometry() const { return _referenceGeometry; }	  

  protected:
    double _W; // Energy Density 

    tvmet::Vector< Vector3D, 3 > _n;  // stress resultants
    tvmet::Vector< Vector3D, 2 > _m;  // moment resultants
    tvmet::Vector< Vector3D, 3 > _nTC;

    std::vector< BaseMaterial_t > _baseMaterials;
  
    ShellGeometry _deformedGeometry;
    ShellGeometry _referenceGeometry;
    double _stretch;

  };

} //namespace voom

#endif //  !defined(__ShellMaterial_h__)
