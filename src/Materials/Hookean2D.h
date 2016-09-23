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
  \file Hookean2D.h

  \brief Concrete class for an isotropic, 2-D (plane stress/strain),
  small-strain, Hookean material object.

*/

#if !defined(__Hookean2D_h__)
#define __Hookean2D_h__

#include "LinearizedMaterial.h"

namespace voom
{

//! Planar material obeying Hooke's law.

/*!  Concrete class for an isotropic, 2-D (plane stress/strain),
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
  indices taking values of 1 and 2.

  The material can obey either a plane strain or plane stress
  formulation depending on the Option passed to the constructor which
  also takes the Young's modulus \f$E\f$ and the Poisson's ratio
  \f$\nu\f$ as arguments.  The Lam&eacute; constants \f$\lambda\f$ and
  \f$\mu\f$ are computed and stored along with \f$E\f$ and \f$\nu\f$
  for convenience (since memory's cheap).

  Stress, strain, and the moduli are represented using Voigt notation
  and stored in arrays (1x3 for stress and strain, 3x3 for moduli) of
  the parent class <tt>LinearizedMaterial<></tt>.
  
*/

class Hookean2D : public LinearizedMaterial<3>
{
public:
  //! Enumeration type to indicate plane stress vs. plane strain behavior.
  enum Option {
    planeStrain,
    planeStress
  }; 

  //! Constructor taking plane strain/stress parameter along with Young's modulus and Poisson's ratio as arguments.
  Hookean2D( Option o, double E, double nu ); 

  //! <tt>LinearizedMaterial<></tt> method implemented to compute the material state.
  void updateState(bool f0, bool f1, bool f2);
  
private:
  //! Young's modulus.
  double _E;
  //! Poisson's ratio.
  double _nu;
  //! First Lam&eacute; constant.
  double _lambda;
  //! Second Lam&eacute; constant (aka, shear modulus).
  double _mu;

  //! Parameter to toggle plane stress vs. plane strain behavior.
  Option _option;
};

} //namespace voom

#endif //  !defined(__Hookean2D_h__)
