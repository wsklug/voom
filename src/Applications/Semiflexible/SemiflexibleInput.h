// Semiflexible Gel Input Class Header

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

   
typedef std::map< std::string, std::string > ParamMap; 


class SemiflexibleInput{

public:

	// Constructor
    SemiflexibleInput(std::string);
    
    // Accessors
	double getReal(std::string);
    int getInt(std::string);
    std::string getStr(std::string);
    
    // Mutators
    void setReal(std::string, double);
    void setInt(std::string, int);
    void setStr(std::string, std::string);
    
    // Destructor 
    // ~SemiflexibleInput();
    
    double qq;

    
private:

	ParamMap pm;
    
};


#endif /* defined(SEMIFLEXIBLEINPUT_H) */
