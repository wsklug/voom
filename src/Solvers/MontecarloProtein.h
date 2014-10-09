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

#if !defined(__MontecarloProtein_h__)
#define __MontecarloProtein_h__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <blitz/array.h>
#include <vector>


#include "ProteinBody.h"
#include "../Applications/Archaea/Utils/PrintingProtein.h"

using namespace std;

namespace voom
{

  class MontecarloProtein
  {
 
  public:

    enum TempSchedule { 
      CONSTANT,         // T = T1
      LINEAR,           // T ~ alpha*t
      FAST,             // T ~ 1/t ("Fast Annealing")
      EXPONENTIAL,      // T ~ alpha^t,  0 < alpha < 1
      STEPWISE
    };
   
    //! Default Constructor
    MontecarloProtein(vector<ProteinNode* > & Proteins,
		      ProteinBody * Body,
		      map<DeformationNode<3> *, vector<DeformationNode<3> *> > & PossibleHosts,
		      int Method,
		      PrintingProtein * PrintProtein,
		      double ResetT = -1.0,
		      int PrintEvery = 1,
		      unsigned int NSteps = 1000,
		      bool print = false): 
      _proteins(Proteins), _body(Body), _possibleHosts(PossibleHosts),
      _method(Method),
      _printProtein(PrintProtein), _resetT(ResetT), _printEvery(PrintEvery),
      _nSteps(NSteps), _print(print), _f(0.0), _fSaved(0.0) {};
    
    //! destructor
    virtual ~MontecarloProtein() {};

    void solve(uint ComputeNeighInterval, double Rsearch);
    bool changeState();
    bool changeAll();

    void SetTempSchedule(TempSchedule Tsched = EXPONENTIAL,
			 const double T01 = 0.0,
			 const double T02 = 0.0,
			 const double FinalTratio = 1.0e-8)
    {
      _Tsched = Tsched;
      _T01 = T01;
      _T02 = T02;
      _FinalTratio = FinalTratio;
    }

    double ComputeUavgSquare(vector<DeformationNode<3>::Point > & OriginalLocations);
    
  private:	

    vector<ProteinNode* > & _proteins;
    ProteinBody * _body;
    map<DeformationNode<3> *, vector<DeformationNode<3> * > > & _possibleHosts;
    int _method;
    PrintingProtein * _printProtein;
    double _resetT;
    int _printEvery;
    unsigned int _nSteps;
    bool _print;

    double _f;
    double _fSaved;

    TempSchedule _Tsched;
    double _T01;
    double _T02;
    double _T1;
    double _T2;
    double _FinalTratio;

  };
  
}; // namespace voom

#endif // __MontecarloProtein_h__
