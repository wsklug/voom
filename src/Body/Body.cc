// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
#include "Body.h"


namespace voom {

  //! check consistency of derivatives
  void Body::checkConsistency() {

    bool verbose=false;

    if(verbose) std::cout << "Body::checkConsistency()" << std::endl;
    blitz::Array<double,1> f_n;
    f_n.resize(_dof);
    f_n = 0.0;

    //
    // analytical residual
    //

    const double h = 1.0e-6;

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();

    for(int a=0, ia=0; a<_nodes.size(); a++) {
      for(int i=0; i<_nodes[a]->dof(); i++, ia++) {
	
	_nodes[a]->addPoint( i, h ); 
	compute(true,false,false);
	f_n(ia) = energy();

	_nodes[a]->addPoint( i, -2.0*h ); 
	compute(true,false,false);
	f_n(ia) -= energy();
	
	_nodes[a]->addPoint( i, h ); 
	
	_nodes[a]->setForce(i,0.0);

      }
    }
    f_n /= 2.0*h;

    compute(false, true, false);

    double Ferror=0.0;
    double Fnorm =0.0;
    double Kerror=0.0;
    double Knorm =0.0;
    double tol=10.0*h;

    std::ofstream ana("f.analytical");
    std::ofstream num("f.numerical");
    std::ofstream dif("f.difference");

    for(int a=0, ia=0; a<_nodes.size(); a++) {
      for(int i=0; i<_nodes[a]->dof(); i++, ia++) {
	double f=_nodes[a]->getForce(i);
	double fn = f_n(ia);
	Ferror = std::max(std::abs(f-fn),Ferror);
	Fnorm += (f)*(f);
	ana << f << " ";
	num << fn << " ";
	dif << f-fn << " ";
      }
      ana << std::endl;
      num << std::endl;
      dif << std::endl;
    }
    Fnorm = sqrt(Fnorm);

    std::cout <<std::setw(18)<<"Fnorm =" <<std::setw(12)<<Fnorm
	      <<std::setw(18)<<"Ferror ="<<std::setw(12)<<Ferror 
	      <<std::setw(18)<<"tol*Fnorm =" <<std::setw(12)<<tol*Fnorm
	      <<std::endl;

    if( Ferror < tol*Fnorm ) {
      std::cout << "Body consistency check PASSED!"<<std::endl;
      return;
    }
    std::cout << "Body consistency check FAILED!"<<std::endl;
    
    return;
    
  }

}
