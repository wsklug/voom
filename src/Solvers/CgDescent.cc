#include "CgDescent.h"
// #include "cg_descent.h"

struct cg_stats {};

extern "C" int cg_descent /*  return  0 (convergence tolerance satisfied)
                           1 (change in func <= feps*|f|)
                           2 (total iterations exceeded maxit)
                           3 (slope always negative in line search)
                           4 (number secant iterations exceed nsecant)
                           5 (search direction not a descent direction)
                           6 (line search fails in initial interval)
                           7 (line search fails during bisection)
                           8 (line search fails during interval update)
                           9 (debugger is on and the function value increases)*/

(
    double      grad_tol , /* StopRule = 1: |g|_infty <= max (grad_tol,
                                            StopFac*initial |g|_infty) [default]
                              StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double            *x , /* input: starting guess, output: the solution */
    int              dim , /* problem dimension (also denoted n) */
    double    (*cg_value)  /* user provided routine to return the function */
               (double *), /* value at x */
    void       (*cg_grad)  /* user provided routine, returns in g the */
      (double *, double*), /* gradient at x*/
    double         *work , /* working array with at least 4n elements */
    double          step , /* initial step for line search
                              ignored unless Step != 0 in cg.parm */
    cg_stats      *Stats   /* structure with statistics(see cg_descent.h) */
) ;


namespace voom {

  int CgDescent::solve(Model * model) {
    ::cgModel = model;
    ::cgSolver = this;

    // resize arrays if necessary
    if(  _x.size() != model->dof() ) {
      resize( model->dof() );
    }
    assert( _x.size() == model->dof() );

    // copy starting guess from model
    model->getField( *this );

    int n = _x.size();
    double step = 1.0e-6;
    int status = 0;
    
    cg_stats Stats;

    Vector_t xtemp;
    xtemp.resize(n);
    xtemp = _x;

    std::cout << "calling cg_descent" << std::endl;
    status = cg_descent( _tol, xtemp.data(), n, 
			 ::myvalue, ::mygrad, 
			 _work.data(), step, &Stats);
    std::cout << "finished cg_descent" << std::endl;

    _x = xtemp;
    model->putField( *this );
    model->computeAndAssemble( *this, true, true, false );

    return 0;
  }


} // end namespace

double myvalue( double *x ) {
  for(int i=0; i < (::cgSolver)->size(); i++) (::cgSolver)->field(i) = x[i];
  ::cgModel->putField( *::cgSolver );
  ::cgModel->computeAndAssemble( *::cgSolver, true, false, false );
  return ::cgSolver->function();
}

void mygrad(double *g, double *x) {
  for(int i=0; i < (::cgSolver)->size(); i++) (::cgSolver)->field(i) = x[i];
  ::cgModel->putField( *::cgSolver );
  ::cgModel->computeAndAssemble( *::cgSolver, false, true, false );
  for(int i=0; i < (::cgSolver)->size(); i++) g[i] = (::cgSolver)->gradient(i);
}
