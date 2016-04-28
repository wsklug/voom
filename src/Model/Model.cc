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
#include<time.h>
#include <fstream>
#include<blitz/array-impl.h>
#include "Model.h"
#include "Solver.h"

// LAPACK FORTRAN subroutine for computing eigenvalues and eigenvectors
extern "C" void dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda,
                       double *w, double *work, int *lwork, int *info);

extern "C" void dsyevx_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *VL, double *VU, int *il, int *iu,
			double *ABSTOL, int *M, double *w, double *Z, int *LDZ, double *work, int *lwork, int *IWORK,
			int *IFAIL, int *info );

extern "C" void dsygvx_( int* ITYPE, char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB,
			 double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
			 double* w, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO );
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
    
    computeAndAssemble(solver,true, true, true);

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

    std::cout<<"Stiffness calculated. Calculating eigenvalues:"<<std::endl;
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
      if ( std::abs(eigenvalues[i]) <= 1.0e-7 )
        computedRank--;
    }
    
    /*for(int i=0; i<_dof; i++) { 
      double h = 0.1/blitz::max(k(i,blitz::Range::all()));
      solver._x += h*k(i,blitz::Range::all());
      putField(solver);
      char name[20];
      sprintf(name,"evec%d",i);
      print(name);
      solver._x -= h*k(i,blitz::Range::all());
    }*/
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
  
  void Model::normalModes(bool eigenvectors, int nmodes, std::string filename) {

    bool verbose = false;
    //for time calculation
    time_t start,end;
    double duration=0.;

    std::cout << "Entered the normal mode calculation." << std::endl;
    Storage solver;
    solver.resize(this->dof());

    blitz::Array<double,2> k_n(_dof,_dof);
    k_n = 0.0;
    getField(solver);
    blitz::Array<double,1>& x = solver._x;

    const double h = 1.0e-8*blitz::max(blitz::abs(x));

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for (int i = 0 ; i < _dof; i++) {
      if (verbose) std::cout << " i = " << i << std::endl;
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

    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;

    // compute Eigenvalues and Eigenvectors by calling LAPACK library
    char jobz;
    if (eigenvectors) jobz = 'V';
    else jobz = 'N';

    char range = 'I';
    char uplo = 'L';
    int n = _dof;
    int lda = n;
    double VL, VU;
    int il=1, iu=nmodes, M;
    double ABSTOL=0.;
    int lwork = 8*n;
    int info;
    double eigenvalues[_dof];
    double work[lwork];
    int IWORK[5*n];
    int IFAIL[n];
    blitz::Array<double,2> Z(nmodes,n);
    int ITYPE=1;
    //stiffness matrix
    blitz::Array<double,2>& k = solver._DDE;

    // calling lapack function here to compute
    // eigenvalues and eigenvectors of inv(m)*k
    
    time(&start);
    dsyevx_( &jobz, &range , &uplo, &n, k.data(), &lda, &VL, &VU, &il, &iu,
     &ABSTOL, &M, eigenvalues, Z.data(), &n, work, &lwork, IWORK,
     IFAIL, &info );
//    dsygvx_(&ITYPE, &jobz, &range, &uplo, &n, k.data(), &lda, mass.data(), &lda,
  //       &VL, &VU, &il, &iu, &ABSTOL, &M,
    //      eigenvalues, Z.data(), &n , work, &lwork, IWORK, IFAIL, &info );
    time(&end);
    duration = difftime(end,start);
    std::cout <<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;
    if (info != 0) {
      std::cout << "Something is wrong in DSYEV_" << std::endl;
      exit(0);
    }
//print out the eigenvectors
    if (eigenvectors) {
      for (int i = 0; i < nmodes; i++) { 
        double h = 1.;//0.3/blitz::max(Z(i,blitz::Range::all()));
        solver._x += h*Z(i,blitz::Range::all());
        putField(solver);
        char name[20];
        std::string fname = filename;
        sprintf(name,"evec%d",i);
        fname += name;
        print(fname);
        solver._x -= h*Z(i,blitz::Range::all());
      }
    }
    std::cout<<"The eigenfrequencies are:"<<std::endl;
    for (int i = 0; i < nmodes; i++) 
      std::cout << "eval[" << std::setw((int)std::ceil(log10(_dof))) << i << "]="
                << eigenvalues[i] << std::endl;
  }
  void Model::normalModes(bool eigenvectors, int nmodes, std::vector<double> &weight, double dens, std::string filename) {

    bool verbose = false;
    //for time calculation
    time_t start,end;
    double duration=0.;

    std::cout << "Entered the normal mode calculation." << std::endl;
    Storage solver;
    solver.resize(this->dof());

    blitz::Array<double,2> k_n(_dof,_dof);
    k_n = 0.0;
    getField(solver);
    blitz::Array<double,1>& x = solver._x;

    const double h = 1.0e-8*blitz::max(blitz::abs(x));

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for (int i = 0 ; i < _dof; i++) {
      if (verbose) std::cout << " i = " << i << std::endl;
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

    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;

    // compute Eigenvalues and Eigenvectors by calling LAPACK library
    char jobz;
    if (eigenvectors) jobz = 'V';
    else jobz = 'N';

    char range = 'I';
    char uplo = 'L';
    int n = _dof;
    int lda = n;
    double VL, VU;
    int il=1, iu=nmodes, M;
    double ABSTOL=0.;
    int lwork = 8*n;
    int info;
    double eigenvalues[_dof];
    double work[lwork];
    int IWORK[5*n];
    int IFAIL[n];
    blitz::Array<double,2> Z(nmodes,n);
    int ITYPE=1;
    //stiffness matrix
    blitz::Array<double,2>& k = solver._DDE;

    //mass matrix
    blitz::Array<double,1> mass(_dof);
    mass=0;
    for(int i=0;i<weight.size();i++)
       {mass(3*i) = mass(3*i+1) = mass(3*i+2) = weight[i]*dens;}
    //transform K -> M^(-1/2)*K*M^(-1/2) where M is lumped mass matrix
    for(int i=0;i<_dof;i++)
       for(int j=0;j<_dof;j++)
          k(i,j)=k(i,j)/sqrt(mass(i))/sqrt(mass(j));

    // calling lapack function here to compute
    // eigenvalues and eigenvectors of inv(m)*k
    
    time(&start);
    dsyevx_( &jobz, &range , &uplo, &n, k.data(), &lda, &VL, &VU, &il, &iu,
     &ABSTOL, &M, eigenvalues, Z.data(), &n, work, &lwork, IWORK,
     IFAIL, &info );
    time(&end);
    duration = difftime(end,start);
    std::cout <<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;
    if (info != 0) {
      std::cout << "Something is wrong in DSYEV_" << std::endl;
      exit(0);
    }
//print out the eigenvectors
    if (eigenvectors) {
      //take off the mass matrix from the eigenvectors
      for (int i=0;i<_dof;i++)
          Z(blitz::Range::all(),i) /= sqrt(mass(i));
      for (int i = 0; i < nmodes; i++) { 
        double h = 1.;//0.3/blitz::max(Z(i,blitz::Range::all()));
        solver._x += h*Z(i,blitz::Range::all());
        putField(solver);
        char name[20];
        std::string fname = filename;
        sprintf(name,"evec%d",i);
        fname += name;
        print(fname);
        solver._x -= h*Z(i,blitz::Range::all());
      }
    }
    std::cout<<"The eigenfrequencies are:"<<std::endl;
    for (int i = 0; i < nmodes; i++) 
      std::cout << "eval[" << std::setw((int)std::ceil(log10(_dof))) << i << "]="
                << eigenvalues[i] << std::endl;
    std::ofstream evals;
    evals.open("evals.dat");
    for (int i = 0; i < nmodes; i++) 
      evals << eigenvalues[i] << std::endl;
    evals.close();
  }
  void Model::normalModes_all(bool eigenvectors, int nmodes, std::vector<double> &weight, double dens,  std::string filename) {

    bool verbose = false;
    //for time calculation
    time_t start, end;
    double duration=0.;

    std::cout << "Entered the normal mode calculation." << std::endl;
    Storage solver;
    std::cout<<"Initializing the stiffness matrix"<<this->dof()<<std::endl;
    solver.resize(this->dof());
    std::cout<<"Initializing the stiffness matrix"<<std::endl;
    blitz::Array<double,2> k_n(_dof,_dof);
    std::cout<<"Stiffness matrix initialized"<<std::endl;
    k_n = 0.0;
    getField(solver);
    blitz::Array<double,1>& x = solver._x;

    const double h = 1.0e-8*blitz::max(blitz::abs(x));

    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for (int i = 0 ; i < _dof; i++) {
      if (verbose) std::cout << " i = " << i << std::endl;
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

    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;

    // compute Eigenvalues and Eigenvectors by calling LAPACK library
    char jobz;
    if (eigenvectors) jobz = 'V';
    else jobz = 'N';

    char uplo = 'L';
    int n = _dof;
    int lda = n;
    int lwork = 3*n-1;
    int info;
    double eigenvalues[_dof];
    double work[lwork];
    time(&start);    
    blitz::Array<double,2>& k = solver._DDE;

    //mass matrix
    blitz::Array<double,1> mass(_dof);
    mass=0;
    for(int i=0;i<weight.size();i++)
       {mass(3*i) = mass(3*i+1) = mass(3*i+2) = weight[i]*dens;}
    //transform K -> M^(-1/2)*K*M^(-1/2) where M is lumped mass matrix
    for(int i=0;i<_dof;i++)
       for(int j=0;j<_dof;j++)
          k(i,j)=k(i,j)/sqrt(mass(i))/sqrt(mass(j));
    //std::cout<<"Mass matrix not used. Only an identity matrix assumed."<<std::endl;
    // calling lapack function here to compute
    dsyev_(&jobz, &uplo, &n, k.data(), &lda, eigenvalues, work, &lwork, &info);
    time(&end);
    duration=difftime(end,start);
    std::cout<<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;

    if (info != 0) {
      std::cout << "Something is wrong in DSYEV_" << std::endl;
      exit(0);
    }
//print out the eigenvectors
    if (eigenvectors) {
      //take off the mass matrix from the eigenvectors
      //for (int i=0;i<_dof;i++)
      //    k(blitz::Range::all(),i) /= sqrt(mass(i));
 
      for (int i = 0; i < nmodes; i++) { 
        double h = 1.0;//0.3/blitz::max(k(i,blitz::Range::all()));
        solver._x += h*k(i,blitz::Range::all());
        putField(solver);
        char name[20];
        std::string fname = filename;
        sprintf(name,"evec%d",i);
        fname += name;
        print(fname);
        solver._x -= h*k(i,blitz::Range::all());
      }
    }
    std::ofstream evals;
    evals.open("evals.dat");
    std::cout<<"The eigenfrequencies are:"<<std::endl;
    for (int i = 0; i < nmodes; i++) 
      evals << eigenvalues[i] << std::endl;
    evals.close();
  }
  void Model::normalModes_all(bool eigenvectors, int nmodes, std::string filename) {
    
    bool verbose = false;
    //for time calculation
    time_t start, end;
    double duration=0.;
    
    std::cout << "Entered the normal mode calculation." << std::endl;
    Storage solver;
    solver.resize(this->dof());
    
    blitz::Array<double,2> k_n(_dof,_dof);
    k_n = 0.0;
    getField(solver);
    blitz::Array<double,1>& x = solver._x;
    
    const double h = 1.0e-8*blitz::max(blitz::abs(x));
    
    // numerical force and stiffness
    blitz::Range all = blitz::Range::all();
    for (int i = 0 ; i < _dof; i++) {
      if (verbose) std::cout << " i = " << i << std::endl;
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
    
    std::cout << "Stiffness matrix calculated. Calculating eigenvalues." << std::endl;
    
    // compute Eigenvalues and Eigenvectors by calling LAPACK library
    char jobz;
    if (eigenvectors) jobz = 'V';
    else jobz = 'N';
    
    char uplo = 'L';
    int n = _dof;
    int lda = n;
    int lwork = 3*n-1;
    int info;
    double eigenvalues[_dof];
    double work[lwork];
    time(&start);    
    blitz::Array<double,2>& k = solver._DDE;
    
    // calling lapack function here to compute
    dsyev_(&jobz, &uplo, &n, k.data(), &lda, eigenvalues, work, &lwork, &info);
    time(&end);
    duration=difftime(end,start);
    std::cout<<"Time taken for eigenvalue decomposition is "<<duration<<" seconds."<<std::endl;
    
    if (info != 0) {
      std::cout << "Something is wrong in DSYEV_" << std::endl;
      exit(0);
    }
    blitz::Array<double,2> a(6,10);
    a=0;
    std::cout<<"Calculating spherical harmonic decomposition of normal modes"<<std::endl;
    blitz::Array<double,1> SHvector(_dof);
    SHvector=0.;
    for(int l=0;l<5;l++){
      for(int m=0;m<2*l+1;m++){
        for(int i=6;i<nmodes;i++){
          //calculate the SH displacement vector
          for(int j=0;j<_nodes.size();j++){
            std::vector<double> dx(3,0.);
            SHdisp(l,m,_nodes[j]->getPoint(0),_nodes[j]->getPoint(1),_nodes[j]->getPoint(2),dx);
            SHvector(3*j)=dx[0]; SHvector(3*j+1)=dx[1]; SHvector(3*j+2)=dx[2];
          }
          //calculate the dot product of eigenvector and the SHvector
          double almi=0.;
          for(int j=0;j<_dof;j++) almi += SHvector(j)*k(i,j);
          a(l,m) += almi/sqrt(eigenvalues[i]); 
        }
      }
    }
    std::cout<<a<<std::endl;
  }
  void Model::SHdisp(int l,int m,double x,double y,double z,std::vector<double>& dx){
    double ur;
    double r=sqrt(x*x+y*y+z*z);
    if(l==0)
      ur=1./2./sqrt(M_PI);
    
    if(l==1){
      if(m==0)
	ur=sqrt(3./4./M_PI)*x/r;
      if(m==1)
	ur=sqrt(3./4./M_PI)*y/r;
      if(m==2)
	ur=sqrt(3./4./M_PI)*z/r;
    }
    
    if(l==2){
      if(m==0)
	ur=1./4.*sqrt(5./M_PI)*(-x*x-y*y+2*z*z)/r/r;
      if(m==1)
	ur=1./2.*sqrt(15./M_PI)*y*z/r/r;
      if(m==2)
	ur=1./2.*sqrt(15./M_PI)*z*x/r/r;
      if(m==3)
	ur=1./2.*sqrt(15./M_PI)*x*y/r/r;
      if(m==4)
	ur=1./4.*sqrt(15./M_PI)*(x*x-y*y)/r/r;
    }
    if(l==3){
      if(m==0)
	ur=1./4.*sqrt(7./M_PI)*z*(2*z*z-3*x*x-3*y*y)/r/r/r;
      if(m==1)
	ur=1./4.*sqrt(35./2./M_PI)*y*(3*x*x-y*y)/r/r/r;
      if(m==2)
	ur=1./4.*sqrt(35./2./M_PI)*x*(x*x-3*y*y)/r/r/r;
      if(m==3)
	ur=1./4.*sqrt(105./M_PI)*z*(x*x-y*y)/r/r/r;
      if(m==4)
	ur=1./2.*sqrt(105./M_PI)*x*y*z/r/r/r;
      if(m==5)
	ur=1./4.*sqrt(21./2./M_PI)*y*(4*z*z-x*x-y*y)/r/r/r;
      if(m==6)
	ur=1./4.*sqrt(21./2./M_PI)*x*(4*z*z-x*x-y*y)/r/r/r;
    }
    
    if(l==4){
      if(m==0)
	ur=3./16.*sqrt(1./M_PI)*(35.*z*z*z*z - 30.*z*z*r*r + 3.*r*r*r*r)/(r*r*r*r);
      if(m==1)
	ur=3./4.*sqrt(5./2./M_PI)*x*z*(7*z*z-3*r*r)/(r*r*r*r);
      if(m==2)
	ur=3./4.*sqrt(5./2./M_PI)*y*z*(7*z*z-3*r*r)/(r*r*r*r);
      if(m==3)
	ur=3./8.*sqrt(5./M_PI)*(x*x-y*y)*(7*z*z-r*r)/(r*r*r*r);
      if(m==4)
	ur=3./4.*sqrt(5./M_PI)*x*y*(7*z*z-r*r)/(r*r*r*r);
      if(m==5)
	ur=3./4.*sqrt(35./2./M_PI)*x*z*(x*x-3*y*y)/(r*r*r*r);
      if(m==6)
	ur=3./4.*sqrt(35./2./M_PI)*x*z*(x*x-3*y*y)/(r*r*r*r);
      if(m==7)
	ur=3./4.*sqrt(35./2./M_PI)*y*z*(3*x*x-y*y)/(r*r*r*r);
      if(m==8)
	ur=3./16.*sqrt(35./M_PI)*(x*x*x*x-6*x*x*y*y+y*y*y*y)/(r*r*r*r);
      if(m==9)
	ur=3./4.*sqrt(35./M_PI)*x*y*(x*x-y*y)/(r*r*r*r);
    }
    
    if(l>4) std::cout<<"Not coded yet"<<std::endl;
    
    
    dx[0]=x/r*ur;
    dx[1]=y/r*ur;
    dx[2]=z/r*ur;
  }
} // namespace voom
