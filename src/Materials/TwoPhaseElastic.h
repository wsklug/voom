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
  \file TwoPhaseElastic.h

  \brief Interface for two phase Evans model for an elastic cytoskeleton & bilayer.

*/

#if !defined(__TwoPhaseElastic_h__)
#define __TwoPhaseElastic_h__

#include<vector>
#include "EvansElastic.h"
#include "VoomMath.h"

 

namespace voom
{

  /*!  Class implementing the two phase
    Evans model for an elastic lipid bilayer & cytoskeleton.
  */
  
  class TwoPhaseElastic : public EvansElastic
  //class EvansElastic : public typename Material< double >
  {
	  
  protected:
    double _kC1;
    double _kC2;
    double _kG1;
    double _kG2;
    double _C01;
    double _C02;
    double _epsilon;//determine the width of the region of transition between phases, length scale
    double _h;      //height of barrier between teh two minima of psi(c)

    double _C; // GL order-parameter field
    Vector2D _dC; // GL order-parameter gradient

    tvmet::Vector< Vector3D, 3 > _nC; // concentration stress resultants
    double _muC; // chemical resultants 
    Vector2D  _muDC;// chemical gradient resultants

  public:

    TwoPhaseElastic(const double kC1, const double kG1, const double C01, const double kC2, const double kG2, const double C02, double mu, const double kS, const double epsilon, const double h) 
      : EvansElastic( kC1, kG1, C01, mu, kS )
    {
      _kC1 = kC1;
      _kG1 = kG1;
      _C01 = C01;
      _kC2 = kC2;
      _kG2 = kG2;
      _C02 = C02;
      _epsilon = epsilon;
      _h = h;
    }

    void updateState(bool f0, bool f1, bool f2);

    void setkCkGC0(){
    
      _kC = _C*_kC1 + (1-_C)*_kC2;
      _kG = _C*_kG1 + (1-_C)*_kG2;
      _C0 = _C*_C01 + (1-_C)*_C02;
    
    }

    void setConcentration(const double& C, const Vector2D & dC) 
    {
      _C = C; 
      _dC = dC;
    }

  };
  
} //namespace voom

#endif //  !defined(__TwoPhaseElastic_h__)
