#if !defined(__MFShape__)
#define __MFShape__

#if (__GNUC__)
#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)
#endif

#include <vector>
#include "voom.h"
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include "Node.h"
#include "NodeBase.h"

// External LAPACK routine to invert matrices 
extern "C" void dgetrf_(int *m, int *n, double *A, int *lda, int *ipiv,
                        int *info);
extern "C" void dgetri_(int *n, double *A, int *lda, int *ipiv, double *work,
                        int *lwork, int *info);

namespace voom
{
  // Class for mesh free shape functions that calculates the shape
  // functions and its derivatives it uses SNNI for calculating the
  // derivatives of shape functions and thus needs a radius
  // corresponding to the volume associated with the point where we need
  // the shape functions and derivatives
  // Three options for kernel - cubic B-spline, singular, modified (primitive+enrichment)
  // May template it to make it suitable for 1D and 2D problems
  
  class MFShape {
  public:
    typedef std::vector<double> FunctionContainer;
    typedef std::vector<int> NodeNContainer;
    typedef tvmet::Vector<double,3> CoordinateArray;
    
    virtual ~MFShape() {}
    
    //temporary constructor
    MFShape(){}
    
    //! return shape functions
    const FunctionContainer & functions() const {return _functions;}
    
    const FunctionContainer & xderivative() const {return _xderivative;}
    
    const FunctionContainer & yderivative() const {return _yderivative;}
    
    const FunctionContainer & zderivative() const {return _zderivative;}
    
    const NodeNContainer & neighb() const {return _nodes;}
    
    //! compute the shape functions and derivatives using cubic B-spline kernel
    void compute(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, const std::vector<int> & list, double radius);
    
    //! compute the shape functions and derivatives using singular kernel
    void compute_sing(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, double radius, double sing_order);
	    
    //! compute the shape functions and derivatives using modified kernel (needs two support sizes)
    void compute_mod(const CoordinateArray & s, const std::vector<DeformationNode<3>* > & nodes, const std::vector<double> & SuppSize, const std::vector<double> & SuppSizeHat, const std::vector<int> & list, double radius);
    
    //! general compute function with argument CoordinateArray instead of DeformationNodes (may delete it later)
    void compute(const CoordinateArray & s, const CoordinateArray* x, const std::vector<double> SuppSize, int nNodes, double radius);
    
    //! cubic B-spline kernel
    double cal_phi(double z){
      //weight function
      /*if(z<1.0){
	if(z>=0.5) return 2*(1-z)*(1-z)*(1-z);
	else return 1-6*z*z+6*z*z*z;
      }
      else return 0;*/
    if(z<1.) return 1.;//-z;
    else return 0.;
    /*  if(z<1.)
        return 2./3.-9./2.*pow(z,2)+19./3.*pow(z,3)-5./2.*pow(z,4);
      else return 0.;*/
    }

    
    //! Singular kernel with 'p' as the order of singularity
    double cal_phi_sing(double z, double p){
      //weight function
    /*  if(z<1.0){
	if(z>=0.5) return 2*(1-z)*(1-z)*(1-z)/pow(z,p);
	else return (1-6*z*z+6*z*z*z)/pow(z,p);
      }
      else return 0;*/
      return cal_phi(z)/pow(z,p);
    }


    bool checkPartitionOfUnity(const CoordinateArray& s, const std::vector<DeformationNode<3>* >& nodes);
    
    //! analytical inverse of a 4*4 matrix
    bool inverse(tvmet::Matrix<double,4,4> M, tvmet::Matrix<double,4,4> & invM);
    
    //! calculates the sum of shape functions at a point using contant basis (for surface defnition)
    float partition_of_unity(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize);
    
    
  protected:
    //! calculates the shape functions at a single point, this is called by compute, returns false if moment matrix is singular, true otherwise
    bool compute_function(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, int nNeigh, FunctionContainer & function);
    
    bool compute_function_sing(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, double sing_order, int nNeigh, FunctionContainer & function);
    
    bool compute_function_mod(const CoordinateArray & s, const std::vector<DeformationNode<3>* >& nodes, const std::vector<double>& SuppSize, const std::vector<double>&  SuppSizeHat, int nNeigh, FunctionContainer & function);
    
    bool compute_function(const CoordinateArray & s, const CoordinateArray* x, const std::vector<double> SuppSize, int nNeigh, FunctionContainer & function);
    
    //
    // //data
    //
    
    //! shape functions (size Node_n)
    FunctionContainer _functions, _xderivative, _yderivative, _zderivative;
    
    //! node number corresponding to the shape functions calculated
    NodeNContainer _nodes;

    // smoothing radius factor
/* Test for GCC < 6.0.0 */
#if (__GNUC__)
#if GCC_VERSION < 60000
    const static double _smooth=0.6;
#else
	static constexpr double _smooth=0.6;
#endif
#else
	static constexpr double _smooth=0.6;
#endif
    //! derivatives of shape functions 
    //DerivativeContainer _derivatives;
  };
}
#endif
