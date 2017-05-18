// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          William S. Klug
//                           Luigi Perotti
//                University of California Los Angeles
//                   (C) 2004 All Rights Reserved
//
//----------------------------------------------------------------------
//
// $Log$
// Revision 1.3  2005/08/22 22:28:28  klug
// Model::setField renamed putField
//
// Revision 1.2  2005/06/27 04:05:11  klug
// *** empty log message ***
//
// Revision 1.1  2005/05/23 18:05:35  klug
// Initial checkin.
//
//----------------------------------------------------------------------

#if !defined(__MontecarloTwoStages_h__)
#define __MontecarloTwoStages_h__

#include "voom.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <blitz/array.h>
#include <vector>
#include "NodeBase.h"
#include "Node.h"
#include "Solver.h"

#include "../Applications/FerroElastic/Utils/PrintingStretches.h"

using namespace std;

namespace voom
{
  /*! Montecarlo solver is applied to a range of dof.
    The rest of the free dof are instead found by energy minimization using CG
  */

  class MontecarloTwoStages: public Solver
  {
 
  public:

    enum TempSchedule { 
      CONSTANT,         // T = T1
      LINEAR,           // T ~ alpha*t
      FAST,             // T ~ 1/t ("Fast Annealing")
      EXPONENTIAL       // T ~ alpha^t,  0 < alpha < 1
    };
   
    //! Default Constructor
    MontecarloTwoStages(vector<vector<ScalarFieldNode<3>* > > & MontecarloDoF,
			vector<int > VarType,
			Solver *GivenSolver,
			PrintingStretches * Printer,
			unsigned int NSteps = 1000,
			bool print = false,
			bool ferroMagnetic = false): 
      _montecarloDoF(MontecarloDoF), _varType(VarType), _solver(GivenSolver), _nSteps(NSteps),  _print(print), _ferroMagnetic(ferroMagnetic), _printingStretches(Printer)
    {
      assert(_montecarloDoF.size() == _varType.size());
      // Compute the number of dof to be changed in the Montecarlo solver
      _size = 0;
      unsigned int i = 0, j = 0;
      for (i = 0; i< _montecarloDoF.size(); i++)
      {
	if (VarType[i] >= 0)
	{
	  _size += _montecarloDoF[i].size();
	}
	  for (j = 0; j< _montecarloDoF[i].size(); j++)
	  {
	    _x.push_back(_montecarloDoF[i][j]->point()); 
	    _xSaved.push_back(_montecarloDoF[i][j]->point());
	  }
      }

      _f = _fSaved = 0.0;
      _flip = 0;

      _printingStretches->setPrintingConfig(true);
    }
    
    //! destructor
    virtual ~MontecarloTwoStages() {};

    //! overloading pure virtual function solve()
    int solve(Model * m);

    void SetTempSchedule(TempSchedule Tsched = FAST,
			 const double T01 = 0.0,
			 const double T02 = 0.0,
			 const double FinalTratio = 1.0e-8)
    {
      _Tsched = Tsched;
      _T01 = T01;
      _T02 = T02;
      _FinalTratio = FinalTratio;
    }
    
    double & field(int i) {return _x[i];}
    double & function() {return _f;}
    double & gradient(int i) {
      std::cerr << "No stiffness in MontecarloTwoStages solver." << std::endl;
      exit(0);
    };    
    double & hessian(int i, int j) {
      std::cerr << "No stiffness in MontecarloTwoStages solver." << std::endl;
      exit(0);
    };
    
    const double field(int i) const {return _x[i];}
    const double function() const {return _f;}
    const double gradient(int i) const {
      std::cerr << "No gradient in MontecarloTwoStages solver." << std::endl;
      exit(0);
    };
    const double hessian(int i, int j) const {
      std::cerr << "No stiffness MontecarloTwoStages solver." << std::endl;
      exit(0);
    };

    double & hessian(int i) {
      std::cerr << "No stiffness MontecarloTwoStages solver." << std::endl;
      exit(0);
    };
    const double hessian(int i) const {
      std::cerr << "No stiffness MontecarloTwoStages solver." << std::endl;
      exit(0);
    };


        
    void zeroOutData(bool f0, bool f1, bool f2) {
      if(f0) _f = 0.0;
    }
    

    void resize(size_t sz) { 
      _x.assign(sz, 0.0); 
      _xSaved.assign(sz, 0.0);
      _size = sz;
      _f = _fSaved = 0.0;
    }
    

    int size() const {return _x.size();}
      
  private:	

    double _f;
    double _fSaved;
    vector<double > _x;
    vector<double > _xSaved;
    
    size_t _size;


    vector<vector<ScalarFieldNode<3>* > > & _montecarloDoF;
    vector<int > _varType;
    Solver * _solver;
    unsigned int _nSteps;
    bool _print;
    bool _ferroMagnetic;

    PrintingStretches *_printingStretches;

    TempSchedule _Tsched;
    double _T01;
    double _T02;
    double _T1;
    double _T2;
    double _FinalTratio;
    int _flip;

    bool changeState(Model *model);

  };
  
}; // namespace voom

#endif // __MontecarloTwoStages_h__
