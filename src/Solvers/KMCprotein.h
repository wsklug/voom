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

#if !defined(__KMCprotein_h__)
#define __KMCprotein_h__

#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>
#include <blitz/array.h>
#include <vector>

// #include <random>
// #include <functional>
// #include <chrono>


#include "ProteinBody.h"
#include "../Applications/Archaea/Utils/PrintingProtein.h"

using namespace std;

namespace voom
{

  class KMCprotein
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
    KMCprotein(vector<ProteinNode* > & Proteins,
		      ProteinBody * Body,
		      map<DeformationNode<3> *, vector<DeformationNode<3> *> > & PossibleHosts,
		      int Method,
		      PrintingProtein * PrintProtein,
		      int PrintEvery = 1,
	              unsigned int NSteps = 1000,
	              int NT = 10,
	              double Zmin = 0.0,
	              double Zmax = 0.0,
		      bool print = false): 
      _proteins(Proteins), _body(Body), _possibleHosts(PossibleHosts),
      _method(Method),
      _printProtein(PrintProtein), _printEvery(PrintEvery),
      _nSteps(NSteps), _NT(NT), _Zmin(Zmin), _Zmax(Zmax), _print(print), 
      _fSaved(0.0), _time(0.0), _approxTime(0.0), _rmax(0.0) 
    {
      _proteinsSize = Proteins.size();
      _rates.assign(_proteinsSize, 0.0);
      // Initialize _temphost vector
      for (int i = 0; i < _proteinsSize; i++) {
	_tempHosts.push_back(_proteins[i]->getHost());
      } 
      /* initialize random seed: */
      int RandSeed = int(time(NULL));
      srand(RandSeed);
      cout << "Seed for random number generator is " << RandSeed << endl;
      // // Setup random number generator
      // _seed = std::chrono::system_clock::now().time_since_epoch().count();
      // _randGenerator(_seed);
      // cout << "Seed for random number generator is " << _seed << endl;
    };
    
    //! destructor
    virtual ~KMCprotein() {};

    void solve(uint ComputeNeighInterval, double Rsearch);
    void computeRates();
    void computeRmax(){ _rmax = *(max_element( _rates.begin(), _rates.end() ) ); };

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

    vector<double > ComputeUavgSquare(vector<DeformationNode<3>::Point > & OriginalLocations);
    
  private:	

    vector<ProteinNode* > & _proteins;
    ProteinBody * _body;
    map<DeformationNode<3> *, vector<DeformationNode<3> * > > & _possibleHosts;
    vector<DeformationNode<3> * > _tempHosts;
    int _method;
    PrintingProtein * _printProtein;
    int _proteinsSize;

    // unsigned int _seed;
    // default_random_engine _randGenerator;

    vector<double > _rates;

    int _printEvery;
    unsigned int _nSteps;
    int _NT;
    double _Zmin;
    double _Zmax;
    bool _print;

    double _fSaved;

    double _time;
    double _approxTime;
    double _rmax;

    TempSchedule _Tsched;
    double _T01;
    double _T02;
    double _T1;
    double _T2;
    double _FinalTratio;

  };
  
}; // namespace voom

#endif // __KMCprotein_h__
