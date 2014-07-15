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
// Revision 1.3  2005/05/23 17:39:11  klug
// Added cvs logging.
//
//----------------------------------------------------------------------

/*! 
  \file StVenant.h

  \brief Interface for a St. Venant-Kirchhoff finite deformation
  elasticity model.

*/

#ifndef _STVENANT_H_
#define _STVENANT_H_

#include "Material.h"
#include "VoomMath.h"

namespace voom {

	class StVenant : public Material
	{

	public:

		//      Constructors/destructors:

		StVenant() { _init(0.0,0.0,0.0); }
		StVenant(double rho, double E, double nu) {
			_init(rho, E, nu);
		}
		virtual ~StVenant() {}
		StVenant(const StVenant &);

		//      Accessors/mutators:

		inline double massDensity();
		inline double longitudinalWaveSpeed();

		//      Operators

		//      General methods:

		void updateState(bool f0, bool f1, bool f2);

		//      Tests:

		void ConsistencyTest();
		static void MFITest();
	        static void IsotropyTest();

	  //void operator=(const StVenant &){;}

	private:

		//      Members:

		double _rho;
		double _E;
		double _nu;
  
	private:
		
	 
		void _init(double rho, double E, double nu);
	};

}
#endif // _STVENANT_H_
