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


    // ***********************************************************************
    // Lin: I cannot believe that you added all of this crap to the
    // shell Material class.  Do all shell materials have
    // concentration, reaction coordinate, etc.?  Obviously not.
    // These things shoulod be only in the material classes that need
    // them.  Don't pollute base classes with unnescessary junk! -WSK
    // ***********************************************************************

    void setConcentration(const double& C, const tvmet::Vector< double, 2 > & dC) {_C = C; _dC = dC;}
    void setReactionCoordinate(const double& eta, const tvmet::Vector< double, 2 > & dEta) {_eta = eta; _dEta = dEta;}

    virtual double energyDensity() const { return _W; }    

    virtual double stretchingEnergy() const = 0;   

    virtual double bendingEnergy() const = 0;

    virtual double concentration() const { return _C; }
  
    virtual const tvmet::Vector< Vector3D, 3 >& stressResultants() const 
    { return _n; }
    virtual const tvmet::Vector< Vector3D, 2 >& momentResultants() const 
    { return _m; }
    virtual const tvmet::Vector< Vector3D, 3 >& concentrationStressResultants() const 
    { return _nC; }
    virtual const double& chemicalResultants() const
    {return _muC;}
    virtual const tvmet::Vector< double, 2 >& chemicalGradientResultants() const
    {return _muDC;}
    virtual const tvmet::Vector< Vector3D, 3 >& totalCurvatureStressResultants() const 
    { return _nTC; }

    virtual const tvmet::Vector< Vector3D, 3 >& GLStressResultants() const 
    { return _nEta; }
    virtual const double& GLResultants() const
    {return _muEta;}
    virtual const tvmet::Vector< double, 2 >& GLGradientResultants() const
    {return _muDEta;}

    const double& stretch() const { return _stretch; }

    const ShellGeometry& shellGeometry() const { return _deformedGeometry; }
    const ShellGeometry& refShellGeometry() const { return _referenceGeometry; }	  

  protected:
    double _W; // Energy Density 

    tvmet::Vector< Vector3D, 3 > _n;  // stress resultants
    tvmet::Vector< Vector3D, 2 > _m;  // moment resultants
    tvmet::Vector< Vector3D, 3 > _nC; // concentration stress resultants
    double _muC;                      // chemical resultants 
    tvmet::Vector< double, 2 >  _muDC;// chemical gradient resultants
    tvmet::Vector< Vector3D, 3 > _nTC;

    tvmet::Vector< Vector3D, 3 > _nEta;
    double _muEta;                      
    tvmet::Vector< double, 2 >  _muDEta;



    std::vector< BaseMaterial_t > _baseMaterials;
  
    ShellGeometry _deformedGeometry;
    ShellGeometry _referenceGeometry;
    double _stretch;

    double _C; //concentration
    tvmet::Vector< double, 2 > _dC;//concentration gradient

    double _eta; //reaction coordinate
    tvmet::Vector< double, 2 > _dEta;//gradient
  };

} //namespace voom

#endif //  !defined(__ShellMaterial_h__)
