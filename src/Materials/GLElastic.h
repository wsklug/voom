// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log: EvansElastic.h,v $
// Revision 1.1  2005/09/28 03:49:23  ma
// Adding cytoskeleton part to SCElasttic, using Evans model
// to claculate cytoskeleton energy
//
// Revision 1.5  2005/05/23 17:38:39  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file GLElastic.h

  \brief Interface for structural phase transformation based on Ginzburg-Landau theory.

*/

#if !defined(__GLElastic_h__)
#define __GLElastic_h__

#include<vector>
#include "FVK.h"
#include "VoomMath.h"

 
namespace voom
{

  /*!  Class implementing the two phase
    Evans model for an elastic lipid bilayer & cytoskeleton.
  */
  
  class GLElastic : public FVK
  //class EvansElastic : public typename Material< double >
  {
    //  public:
    //typedef tvmet::Vector<double,3> Vector3D;
    //typedef tvmet::Matrix<double,2,2> Tensor2D;
    //typedef tvmet::Matrix<double,3,3> Tensor3D;
	  
  public:

    GLElastic( const double kC, 
	       const double kG, 
	       const double C0, 
	       const double E, 
	       const double nu, 
	       const double Gamma,
	       const double g0,
	       const double deltag,
	       const double muA,
	       int formulation=0) 
      : FVK(kC, kG, C0, E, nu), 
	_Gamma(Gamma), _g0(g0), _deltag(deltag), _muA(muA),
	_formulation(formulation)
    { 
      assert(_formulation==0 || _formulation==1);
    }

    void updateState(bool f0, bool f1, bool f2);
    
    virtual double Field() const { return _eta; }
    
    void setField(const double eta, const Vector2D & deta) 
    {_eta = eta; _deta = deta;}

    virtual double GLResultants() const {return _mu;}
    virtual const Vector2D & GLGradientResultants() const {return _lambda;}

  protected:
    double _Gamma;//determine the width of the region of transition between phases, length scale
    double _g0;      //curvature of the reaction coordinate free energy
    double _deltag; // difference in energy eta 0->1

    double _muA; // Constant chemical potential
    double _eta; // GL order-parameter field
    Vector2D _deta; // GL order-parameter gradient

    double _mu; // chemical resultants 
    Vector2D _lambda; // chemical gradient resultants

    int _formulation; // what type of GL function is it?
  };
  
} //namespace voom

#endif //  !defined(__GLElastic_h__)
