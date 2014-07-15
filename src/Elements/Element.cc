// -*- C++ -*-
//----------------------------------------------------------------------
//
//                         William S. Klug
//                University of California Los Angeles
//                 (C) 2004-2005 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.2  2005/10/21 01:07:21  klug
// Cleaned and removed old stuff related to storing forces and stiffness
// in elements.  Note: stiffness part needs to be fixed in checkConsistency()
// and checkRank() is totally broken.
//
// Revision 1.1  2005/08/22 22:18:30  klug
// Assembly shifted from elements to nodes.  Element class renamed
// StandardElement, ElementBase class renamed Element, and verification
// tests moved to Element class and their definitions moved from Element.icc
// to Element.cc which is now compiled into libElement.a.
//
// Revision 1.2  2005/05/23 17:25:08  klug
// Added index map; renamed strainEnergy as energy.
//
// Revision 1.1  2005/04/11 16:24:12  klug
// New Element base classes with forces and stiffness stored locally.
//
//----------------------------------------------------------------------

/*! 
  \file Element.cc

  \brief Class for a finite element.

*/

#include "Element.h"
#include "VoomMath.h"

// LAPACK FORTRAN subroutine for computing eigenvalues and eigenvectors
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                       double *w, double *work, int *lwork, int *info);

namespace voom {

  bool Element::checkConsistency() {

    srand(time(0));

    int a=0, ai=0;
    for( BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++,a++ ){
      for( int i=0; i<(*n)->dof(); i++, ai++) {
	(*n)->setForce(i,0.0);
      }
    }

    std::cout << "Total number of DOF: " << ai << std::endl;

    blitz::Array<double,1> forces_n(ai);

    double h = 1.0e-8;
    a=ai=0;
    for( BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++,a++ ){
      for( int i=0; i<(*n)->dof(); i++, ai++) {
        // perturb +
	(*n)->addPoint(i,h);

        compute(true,false,false);

        forces_n(ai) = _energy;

        // (*n)->work() ;
                                
        // perturb -
	(*n)->addPoint(i,-h-h);

        compute(true,false,false);

        forces_n(ai) -= _energy;

        //      + (*n)->work() ;
        forces_n(ai) /= 2.0*h;

	(*n)->addPoint(i,h);

      }
    }
   
    compute(false,true,false);
    double Ferror=0.0;
    double Fnorm =0.0;
    double tol=100.0*h;
    a=ai=0;
    for( BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++, a++ ) {
      for( int i=0; i<(*n)->dof(); i++, ai++) {
	const double f = (*n)->getForce(i);
	const double fn = forces_n(ai);
	Ferror = std::max(std::abs(f-fn),Ferror);
	Fnorm += (f)*(f);
      }
    }
    Fnorm = sqrt(Fnorm);
    
    std::cout << std::setw(8) << "dof" 
	      << std::setw(16) << "f"
	      << std::setw(16) << "fn"
	      << std::endl;
    a=ai=0;
    for( BaseNodeIterator n=_baseNodes.begin(); n!=_baseNodes.end(); n++, a++ ) {
      for( int i=0; i<(*n)->dof(); i++, ai++) {
	const double f = (*n)->getForce(i);
	const double fn = forces_n(ai);
	std::cout << std::setw(8) << ai
		  << std::setw(16) << f
		  << std::setw(16) << fn
		  << std::endl;
      }
    }
    std::cout <<std::setw(18)<<"Ferror ="<<std::setw(12)<<Ferror 
	      <<std::setw(18)<<"Fnorm ="<<std::setw(12)<<Fnorm
	      <<std::setw(18)<<"tol*Fnorm =" <<std::setw(12)<<tol*Fnorm
	      <<std::endl;
    if( Ferror < tol*Fnorm ) {
      std::cout << "Element consistency check PASSED!"<<std::endl;
      return true;
    }
    std::cout << "Element consistency check FAILED!"<<std::endl;

    return false;
  }

  bool Element::checkRank(const int rank) {
 //    compute(false,false,false,true);
//     srand(time(0));
//     const int nDOF = _field.size();
//     for(int ai=0; ai<nDOF; ai++) 
//       _field(ai) += _field(ai)*0.05*static_cast<double>(rand())/RAND_MAX;

//     compute(false,false,true);
    
//     // compute Eigenvalues and Eigenvectors by calling LAPACK library
//     char jobz = 'V';
//     char uplo = 'L';
//     int  n    = nDOF;
//     int  lda  = n;
//     int  lwork = 3*n-1;
//     int  info;
//     double eigenvalues[nDOF];
//     double work[lwork];
    
//     ElementMatrix k(_stiffness.shape());
//     k = _stiffness;
//     // calling lapack function here to compute
//     // eigenvalues and eigenvectors of k
//     dsyev_(&jobz, &uplo, &n, k.data(),
//            &lda, eigenvalues, work, &lwork, &info);


//     if (info != 0) {
//       std::cout << "Something is wrong in DSYEV_" << std::endl;
//       exit(0);
//     }

//     int computedRank = nDOF;
//     for (int i = 0; i < nDOF; i++){
//       if ( std::abs(eigenvalues[i]) <= 1.0e-8 )
//         computedRank--;
//     }

//     std::cout << "Rank as computed: " << computedRank << std::endl
// 	      << "as it should be:  " << rank << std::endl;

//     if( computedRank != rank ) {
//       std::cout << "Element rank test FAILED!!!!"<<std::endl;
//       return false;
//     }

//     std::cout << "Element rank test PASSED!!!!"<<std::endl;
//     return true;
    return false;
  }

}; // namespace voom
