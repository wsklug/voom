#if !defined(__LMEtetShape__)
#define __LMEtetShape__

#include <vector>
#include "voom.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Node.h"
#include "NodeBase.h"
#include "VoomMath.h"

using namespace std;

// External LAPACK routine to invert matrices 
extern "C" void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
extern "C" void dpotrs_(char *uplo, int *n, int *nrhs, double *a, int *lda,
                        double *b, int *ldb, int *info);
extern "C" void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv,
                        double *b, int *ldb, int *info);
extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    // multiply two matrices
    void dgemm_(char *transa, char *transb, int* m, int* n, int* k,
		double* alpha, double* a, int* lda, 
		double* b, int* ldb, 
		double* beta, double* c, int* ldc);
}

namespace voom
{
  // Class for Linear Maximum Entropy shape functions (LME)
  // One shape function instantion per each tet element in the pseudo mesh of the domain
  // In this case the shape function acts more like and element
  
  class LMEtet {
  public:
    typedef std::vector<double> FunctionContainer;
    typedef std::vector<int> NodeNContainer;
    typedef tvmet::Vector<double,3> CoordinateArray;
    
    virtual ~LMEtet() {}
    
    //Temporary constructor
    LMEtet(double beta, double searchR, double tol, int nItMax, vector<int > NodesInd, vector<tvmet::Vector<double,3 > > EvaluationPoints, vector<double > Weights);
    
    //! Return shape functions and their derivatives
    const FunctionContainer & functions() const {return _functions;}
    
    const FunctionContainer & xderivative() const {return _xderivative;}
    
    const FunctionContainer & yderivative() const {return _yderivative;}
    
    const FunctionContainer & zderivative() const {return _zderivative;}
    
    const NodeNContainer & neighb() const {return _nodesInd;}
    NodeNContainer & neighb() {return _nodesInd;}


    //! Compute the shape functions and derivatives
    bool compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& Nodes);

    double checkPartitionOfUnity();

    vector<double>  LME_r(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes);
    double LME_pa(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, int);
    vector<vector<double > > LME_J(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes);
   

    vector<vector<double > > DyadicProduct(const vector<double > & a, const vector<double > & b);
    vector<vector<double > > ScalarProduct(const double a, const vector<vector<double > > A);
    void MatrixSum (vector<vector<double > > & A, const vector<vector<double > > & B);
    void MatrixSub (vector<vector<double > > & A, const vector<vector<double > > & B);

    void inverse(double* A, int N)
    {
      int *IPIV = new int[N+1];
      int LWORK = N*N;
      double *WORK = new double[LWORK];
      int INFO;
      
      dgetrf_(&N,&N,A,&N,IPIV,&INFO);
      dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
      
      delete IPIV;
      delete WORK;
    }

    void multiplyAKAT(double* A, int N, double* Kbar, double* K)
    {
      double *WORK = new double[N*N];
      char trN = 'N';
      char trY = 'T';
      double alpha = 1.0, beta = 0.0;
      
      dgemm_(&trN, &trN, &N, &N, &N,
	     &alpha, A, &N,
	     Kbar, &N,
	     &beta, WORK, &N);
      
      dgemm_(&trN, &trY, &N, &N, &N,
	     &alpha, WORK, &N,
	     A, &N,
	     &beta, K, &N);
      
      delete WORK;
    }


    void LocalStiffness(const std::vector<DeformationNode<3>* >& Nodes,  const double Kmat[3][3][3][3], double* LMEK);



  protected:
    // Data
    
    //! Shape functions and their derivatives
    FunctionContainer _functions, _xderivative, _yderivative, _zderivative;
    
    //! Node number corresponding to the shape functions calculated
    NodeNContainer _nodesInd;

    // Quadrature points position in the tet element
    vector<tvmet::Vector<double,3 > > _evaluationPoints;
    vector<double > _weights;

    // Search radius, parameter beta, minimizer lambda
    double _beta;
    double _searchR;
    double _tol;
    int _nItMax;
    vector<double> _lambda;

  };
}
#endif
