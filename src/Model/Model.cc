// -*- C++ -*-
//----------------------------------------------------------------------
//
//                    William S. Klug & Feng Feng
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.13  2005/08/22 22:30:52  klug
// Assembly shifted from elements to nodes.  Bodies now compute energy.
// Model::setField renamed less ambiguously to putField.
//
// Revision 1.12  2005/05/23 17:54:15  klug
// New interface with data storage pushed to solver.
//
//----------------------------------------------------------------------

#include <fstream>
#include<blitz/array-impl.h>
#include "Model.h"
#include "Solver.h"

// LAPACK FORTRAN subroutine for computing eigenvalues and eigenvectors
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                       double *w, double *work, int *lwork, int *info);


namespace voom
{

  Model::Model( const BodyContainer & bodies, const NodeContainer & nodes )
  {
    _bodies = bodies;
    _nodes = nodes;
    std::cout << std::setw(15)<<"Building Model from "
	      <<_bodies.size()
	      << " Bodies and "
	      <<_nodes.size()
	      << " Nodes." << std::endl;
    _dof = 0;
    // The model used to get nodes from the bodies.  Now the master
    // list of nodes is stored in the model, so we count dof directly.
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
      _dof += (*n)->dof();
    std::cout << "Created Model with "<< _dof << " dof."<<std::endl;

#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    
  }

  Model::Model( const NodeContainer & nodes )
  {
    _nodes = nodes;
    std::cout << std::setw(15)<<"Building Model from "
	      <<_nodes.size()
	      << " Nodes." << std::endl;
    _dof = 0;
    // The model used to get nodes from the bodies.  Now the master
    // list of nodes is stored in the model, so we count dof directly.
    for(ConstNodeIterator n=_nodes.begin(); n!=_nodes.end(); n++) 
      _dof += (*n)->dof();
    std::cout << "Created Model with "<< _dof << " dof."<<std::endl;

#ifdef WITH_MPI
    MPI_Comm_size( MPI_COMM_WORLD, &_nProcessors );
    MPI_Comm_rank( MPI_COMM_WORLD, &_processorRank );
#endif
    
  }

  //! check consistency of derivatives
  bool Model::checkConsistency(bool f1, bool f2) {

    bool verbose=false;
    //bool verbose=true;

    std::cout << "Model: checking consistency" << std::endl;

    if(verbose) std::cout << "Model::checkConsistency()" << std::endl;
    Storage solver;
    solver.resize(this->dof());
    blitz::Array<double,1> f_n;
    if(f1) {
      f_n.resize(_dof);
      f_n = 0.0;
    }
    blitz::Array<double,2> k_n;
    if(f2) {
      k_n.resize(_dof,_dof);
      k_n = 0.0;
    }

    getField(solver);
    blitz::Array<double,1> & x = solver._x;

    if(verbose) std::cout << "solver.field() = " << x << std::endl;

    //
    // analytical residual
    //

    const double h = 1.0e-6*blitz::max(blitz::abs(x));

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for( int i = 0 ; i < _dof; i++) {

      if(verbose) std::cout << "i = " << i;
      else if(i%1000 == 0) std::cout << "i = " << i << std::endl;
      x(i) += h; 
      putField(solver);
      computeAndAssemble(solver, f1, f2, false);
      if(verbose) std::cout << " E+ = " << solver._E;
      if(f1) f_n(i)= solver._E;
      if(f2) k_n(i,all) = solver._DE;

      x(i) -= 2.0*h; putField(solver);
      computeAndAssemble(solver,f1, f2, false);
      if(verbose) std::cout << " E- = " << solver._E;
      if(f1) f_n(i) -= solver._E;
      if(f2) k_n(i,all) -= solver._DE;
      if(verbose) std::cout << " dE = " << f_n(i) << std::endl;
 
      x(i) += h; putField(solver);
    }
    if(f1) f_n /= 2.0*h;
    if(f2) k_n /= 2.0*h;

    //computeAndAssemble(solver,true, true, true);
    computeAndAssemble(solver,true,f1,f2);

    double Ferror=0.0;
    double Fnorm =0.0;
    double Kerror=0.0;
    double Knorm =0.0;
    double tol=10.0*h;
    for(int i=0; i<_dof; i++){
      //std::cout << "i = " << i << "; ";
      const double & f = solver._DE(i);
      //std::cout << "f = " << f << "; ";
      const double & fn = f_n(i);
      //std::cout << "fn = " << fn << std::endl;
      if(f1) {
	Ferror = std::max(std::abs(f-fn),Ferror);
	Fnorm += (f)*(f);
      }
      if(f2) {
	for(int j=0; j<_dof; j++){
	  const double & k = solver._DDE(i,j);
	  const double & kn = k_n(i,j);
	  Kerror = std::max(std::abs(k-kn),Kerror);
	  Knorm += (k)*(k);
	}
      }
    }
    if(f1) Fnorm = sqrt(Fnorm);
    if(f2) Knorm = sqrt(Knorm);
    
    if(f1)
      std::cout <<std::setw(18)<<"Fnorm =" <<std::setw(12)<<Fnorm
		<<std::setw(18)<<"Ferror ="<<std::setw(12)<<Ferror 
		<<std::setw(18)<<"tol*Fnorm =" <<std::setw(12)<<tol*Fnorm
		<<std::endl;
    if(f2)
      std::cout <<std::setw(18)<<"Knorm =" <<std::setw(12)<<Knorm
		<<std::setw(18)<<"Kerror ="<<std::setw(12)<<Kerror 
		<<std::setw(18)<<"tol*Knorm =" <<std::setw(12)<<tol*Knorm
		<<std::endl;
#ifdef WITH_MPI
    if(f1) {
      char name[20];
      sprintf(name,"f.%d",_processorRank);
	  std::ofstream fout(name);
	  fout << solver._DE << std::endl
	       << f_n << std::endl;
    }
#endif
    if( (!f1 || Ferror < tol*Fnorm) && (!f2 || Kerror < tol*Knorm) ) {
      std::cout << "Model consistency check PASSED!"<<std::endl;
      return true;
    }
    std::cout << "Model consistency check FAILED!"<<std::endl;
    if(f1) {
#ifdef WITH_MPI
      if(_processorRank==0)
#endif
	{
	  if(verbose) {
	    std::cout <<"forces = "<<solver._DE<<std::endl
		      <<"forces_n = "<<f_n<<std::endl;
	  }
	  std::ofstream ana("f.analytical");
	  ana << solver._DE;
	  std::ofstream num("f.numerical");
	  num << f_n;
	  std::ofstream dif("f.difference");
	  f_n -= solver._DE;
	  dif << f_n;
	}
    }
    if(f2) {
#ifdef WITH_MPI
      if(_processorRank==0)
#endif
	{
	  if(verbose) {
	    std::cout <<"stiffness = "<<solver._DDE<<std::endl
		      <<"stiffness_n = "<<k_n<<std::endl;
	  }
	}
    }
    return false;

  }

  bool Model::checkRank(const int rank, bool numerical) {

    Storage solver;
    solver.resize(this->dof());
    if(numerical) {

      blitz::Array<double,2> k_n(_dof,_dof);

      k_n = 0.0;

      getField(solver);
      blitz::Array<double,1> & x = solver._x;

      //std::cout << x << std::endl;

      //
      // analytical residual
      //

      const double h = 1.0e-8*blitz::max(blitz::abs(x));

      // numerical force and stiffness
      blitz::Range all = blitz::Range::all();
      for( int i = 0 ; i < _dof; i++) {

	x(i) += h; 
	putField(solver);
	computeAndAssemble(solver, false, true, false);
	k_n(i,all) = solver._DE;

	x(i) -= 2.0*h; putField(solver);
	computeAndAssemble(solver, false, true, false);
	k_n(i,all) -= solver._DE;
  
	x(i) += h; putField(solver);
      }
      k_n /= 2.0*h;

      solver._DDE = 0.5*(k_n + k_n.transpose(blitz::secondDim,blitz::firstDim));

    } else {
      computeAndAssemble(solver,false,false,true);
    }
    // compute Eigenvalues and Eigenvectors by calling LAPACK library
    char jobz = 'V';
    char uplo = 'L';
    int  n    = _dof;
    int  lda  = n;
    int  lwork = 3*n-1;
    int  info;
    double eigenvalues[_dof];
    double work[lwork];
    
    blitz::Array<double,2> & k = solver._DDE;
    // calling lapack function here to compute
    // eigenvalues and eigenvectors of k
    dsyev_(&jobz, &uplo, &n, k.data(),
           &lda, eigenvalues, work, &lwork, &info);


    if (info != 0) {
      std::cout << "Something is wrong in DSYEV_" << std::endl;
      exit(0);
    }

    int computedRank = _dof;
    for (int i = 0; i < _dof; i++){
      if ( std::abs(eigenvalues[i]) <= 1.0e-8 )
        computedRank--;
    }
    
    for(int i=0; i<_dof; i++) { 
      double h = 0.1/blitz::max(k(i,blitz::Range::all()));
      solver._x += h*k(i,blitz::Range::all());
      putField(solver);
      char name[20];
      sprintf(name,"evec%d",i);
      print(name);
      solver._x -= h*k(i,blitz::Range::all());
    }
    for(int i=0; i<_dof; i++) 
      std::cout << "eval["<<std::setw((int)std::ceil(log10(_dof)))<<i<<"]="
		<<eigenvalues[i]
		<<std::endl;
    
    std::cout << "Rank as computed: " << computedRank << std::endl
	      << "as it should be:  " << rank << std::endl;

    if( computedRank != rank ) {
      std::cout << "Element rank test FAILED!!!!"<<std::endl;
     return false;
    }

    std::cout << "Element rank test PASSED!!!!"<<std::endl;
    return true;

  }
  
} // namespace voom


