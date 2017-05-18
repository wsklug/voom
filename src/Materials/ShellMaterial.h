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

    virtual const Tensor3D & cauchyStress() { return _cauchy; }	 
    virtual const std::vector<double > invariants() { return std::vector<double >(2, 0.0); }	 

  protected:
    double _W; // Energy Density 

    tvmet::Vector< Vector3D, 3 > _n;  // stress resultants
    tvmet::Vector< Vector3D, 2 > _m;  // moment resultants
    tvmet::Vector< Vector3D, 3 > _nTC;

    std::vector< BaseMaterial_t > _baseMaterials;
  
    ShellGeometry _deformedGeometry;
    ShellGeometry _referenceGeometry;
    double _stretch;

    Tensor3D _cauchy; // Cauchy Stress Tensor
  };

} //namespace voom

#endif //  !defined(__ShellMaterial_h__)
