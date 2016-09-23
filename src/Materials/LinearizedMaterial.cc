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
// Revision 1.5  2005/10/22 19:16:14  klug
// GCC 4 compatibility fix.
//
// Revision 1.4  2005/05/23 17:37:31  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

#include "LinearizedMaterial.h"
#include "VoomMath.h"

namespace voom {

  // Isotropy check for 2-D elasticity
  template<>
  bool LinearizedMaterial<3>::checkIsotropy() {
    // Generate a random rotation
    srand(time(0));
    double theta = 2.0*M_PI*rand()/static_cast<double>(INT_MAX);
    tvmet::Matrix<double,2,2> R;
    R = cos(theta), -sin(theta),
        sin(theta),  cos(theta); 

    // Create random strain, compute w, sigma, c
    VoigtVector epsilon;
    for(int i=0; i<3; i++) {
      epsilon(i) = 0.2 * ( rand()/static_cast<double>(RAND_MAX) - 0.5 );
    }
    setStrain(epsilon);
    updateState(true,true,true);
    double w = energyDensity();
    VoigtVector sigma;
    sigma = stress();
    VoigtMatrix c;
    c = moduli();

    // Rotate strain, compare rotated w, sigma, and c to originals
    tvmet::Matrix<double,2,2> e;
    e =     epsilon(0), 0.5*epsilon(2),
        0.5*epsilon(2),     epsilon(1);
    tvmet::Matrix<double,2,2> e_rot;
    e_rot = tvmet::trans(R)*e*R;
    VoigtVector epsilon_rot;
    epsilon_rot = e_rot(0,0), e_rot(1,1), 2.0*e_rot(0,1);
    setStrain(epsilon_rot);
    updateState(true,true,true);

    double w_rot = energyDensity();
    VoigtVector sigma_rot1;
    sigma_rot1 = stress();
    VoigtMatrix c_rot;
    c_rot = moduli();
    
    tvmet::Matrix<double,2,2> s_rot1;
    s_rot1 = sigma_rot1(0), sigma_rot1(2),
             sigma_rot1(2), sigma_rot1(1);

    // compare
    double wError = fabs((w_rot-w)/w);

    tvmet::Matrix<double,2,2> s;
    s = sigma(0), sigma(2),
        sigma(2), sigma(1);
    tvmet::Matrix<double,2,2> s_rot2;
    s_rot2 = tvmet::trans(R)*s*R;
    VoigtVector sigma_rot2;
    sigma_rot2 = s_rot2(0,0), s_rot2(1,1), s_rot2(0,1);
    double sigmaError = 
       norm2(sigma_rot2 - sigma_rot1)/ norm2(sigma_rot2);

    double tol = 1.0e-12;
    if( wError > tol  || sigmaError > tol ) {
      std::cout<<"LinearizedMaterial::checkIsotropy() FAILED!" 
	       << std::endl;    
      std::cout <<"wError = "<<wError<<std::endl
		<<"sigmaError = "<<sigmaError<<std::endl
		<<"e =" << std::endl
		<< e << std::endl
		<<"e_rot =" << std::endl
		<< e_rot << std::endl
		<<"R =" << std::endl
		<< R << std::endl
		 <<"w =" << std::endl
		<< w << std::endl
		<<"w_rot =" << std::endl
		<< w_rot << std::endl
		<<"s =" << std::endl
		<< s << std::endl
		<<"s_rot1 =" << std::endl
		<< s_rot1 << std::endl
		<<"s_rot2 =" << std::endl
		<< s_rot2 << std::endl;
      return false;
    }
    std::cout<<"LinearizedMaterial::checkIsotropy() PASSED!" 
	     << std::endl;    
    return true;
  }

};
