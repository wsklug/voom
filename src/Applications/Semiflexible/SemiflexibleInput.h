// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          
//                
//                   
//
//----------------------------------------------------------------------

#ifndef SEMIFLEXIBLEINPUT_H
#define SEMIFLEXIBLEINPUT_H

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
//#include <complex>
#include <map> // ?
#include <cmath> // maybe not needed
#include <sstream>
#include <iomanip> //?
#include <blitz/array.h> //?
#include "voom.h"
#include "SemiflexibleGel.h"


#include "GelOutput.h"
#include <algorithm>
#include <cstdio>
#include <ctime>
#include "Node.h"
#include "Lbfgsb.h"


class SemiflexibleInput{
  
 public:

  // Constructor
  SemiflexibleInput(std::string paramFileName);
  
  // Accessors
  double getReal(std::string name) const;
  int getInt(std::string name) const;
  std::string getStr(std::string name) const;
  
  //SemiflexibleGel<2> * gel() const {return _gel;}

  //PeriodicBox * box() const {return _box;}

  // Mutators
  void setReal(std::string name, double value);
  void setInt(std::string name, int value);
  void setStr(std::string name, std::string value);
  
  // Destructor 
  // ~SemiflexibleInput();

  
  
 private:
  
  std::map< std::string, std::string > _pm;

  //SemiflexibleGel<2> * _gel;

  //PeriodicBox * _box;
    
};


#endif /* defined(SEMIFLEXIBLEINPUT_H) */
