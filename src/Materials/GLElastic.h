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

 
typedef tvmet::Vector<double,3> Vector3D;
typedef tvmet::Matrix<double,2,2> Tensor2D;
typedef tvmet::Matrix<double,3,3> Tensor3D;

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
	  
  protected:
    double _Gamma;//determine the width of the region of transition between phases, length scale
    double _gDoublePrime;      //curvature of the reaction coordinate free energy
    double _gammaZero; //coupling constant
    double _factor; //determines the magnitude of two minima, _factor=1 gives two equal


  public:

    GLElastic(const double kC, const double kG, const double C0, const double E, const double nu, 
	      const double gDoublePrime, const double Gamma, const double gammaZero, const double factor) 
      : FVK(kC, kG, C0, E, nu)
    {
      _Gamma = Gamma;
      _gDoublePrime = gDoublePrime;
      _gammaZero = gammaZero;
      _factor = factor;
    }

      void updateState(bool f0, bool f1, bool f2);


  };
  
} //namespace voom

#endif //  !defined(__GLElastic_h__)
