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
// Revision 1.3  2005/05/23 17:36:54  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file Hookean.h

  \brief Concrete class for an isotropic, 3-D, small-strain, Hookean
  material object.

*/

#if !defined(__Hookean_h__)
#define __Hookean_h__

#include "LinearizedMaterial.h"

namespace voom
{

//! 3-D material obeying Hooke's law.

/*!  Concrete class for an isotropic, 3-D,
  small-strain, material object obeying Hooke's law, i.e., 
  \f[
  \sigma_{ij} = c_{ijkl}\epsilon_{kl} 
  \f] 
  where \f$\epsilon_{kl}\f$ are the components of the small
  (infinitesimal) strain tensor, \f$\sigma_{ij}\f$ are the components 
  of the Cauchy stress tensor, and \f$c_{ijkl}\f$ are the elastic moduli, 
  given in terms of the Lam&eacute; constants \f$\lambda\f$ and
  \f$\mu\f$ by
  \f[
  c_{ijkl} = \lambda\delta_{ij}\delta_{kl} + 
  \mu(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk} )
  \f]
  (due to isotropy).
  Summation according to the Einstein convention is assumed with
  indices taking values of 1 and 3.

  The only available constructor takes the Young's modulus \f$E\f$ and
  the Poisson's ratio \f$\nu\f$ as arguments.  The Lam&eacute;
  constants \f$\lambda\f$ and \f$\mu\f$ are computed and stored along
  with \f$E\f$ and \f$\nu\f$ for convenience (since memory's cheap).

  Stress, strain, and the moduli are represented using Voigt notation
  and stored in arrays (1x6 for stress and strain, 6x6 for moduli) of
  the parent class <tt>LinearizedMaterial<></tt>.
*/

class Hookean : public LinearizedMaterial<6>
{
public:
  //! Constructor taking Young's modulus and Poisson's ratio as arguments.
  Hookean( double E, double nu ); 

  //! <tt>LinearizedMaterial<></tt> method implemented to compute the material state.
  virtual void updateState(bool f0, bool f1, bool f2);
  
private:
   //! Young's modulus.
  double _E;
  //! Poisson's ratio.
  double _nu;
  //! First Lam&eacute; constant.
  double _lambda;
  //! Second Lam&eacute; constant (aka, shear modulus).
  double _mu;
};

} //namespace voom

#endif //  !defined(__Hookean_h__)
