//--------------------------------------------------------------------
//
//                  William Klug and Michael Ortiz
//                California Institute of Technology
//                  (C) 2002 - All Rights Reserved
//
//--------------------------------------------------------------------
//
// Rosenbrock.h:  Interface to the Rosenbrock class.
//
//--------------------------------------------------------------------

#ifndef _CDNA_ROSENBROCK_H_
#define _CDNA_ROSENBROCK_H_

#include "Body.h"

namespace voom 
{

  class Rosenbrock : public Body  {

  public:
    // ---------------
    // Public methods:
    // ---------------  

    // 
    // Constructors
    //
    Rosenbrock( const int dim, const bool debug );

    //! destructor
    virtual ~Rosenbrock() {}

    //
    // Output current state of the model
    //
    virtual void printState();

    virtual void compute( bool f0, bool f1, bool f2 );

    virtual double getEnergy() {return _f;}

    //! Copy positions
    virtual void getPositions( blitz::Array< double, 1 > & x ) const;
    //! Copy velocities
    virtual void getVelocities( blitz::Array< double, 1 > & v );
    //! Copy accelerations
    virtual void getAccelerations( blitz::Array< double, 1 > & v );
    //! Copy residual forces
    virtual void getResidual( blitz::Array< double, 1 > & r );
    //! Copy tangent stiffness
    virtual void getStiffness( blitz::Array< double, 1 > & k, 
			       blitz::Array< int, 2 > & ndx );

    // Copy field values from arrays into nodes

    //! Copy positions
    virtual void setPositions( const blitz::Array< double, 1 > & x );
    //! Copy velocities
    virtual void setVelocities( const blitz::Array< double, 1 > & v );
    //! Copy accelerations
    virtual void setAccelerations( const blitz::Array< double, 1 > & a );


  private:
    // ---------------------------------------------
    // Private member data variables and containers:
    // ---------------------------------------------  
    double _f;
    blitz::Array< double, 1 > _x;
    blitz::Array< double, 1 > _grad;
    double _alpha;

    bool _debug;

  }; // class Rosenbrock  

}; // namespace voom

#endif // _ROSENBROCK_H_
