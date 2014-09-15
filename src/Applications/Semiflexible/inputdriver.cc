// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          
//                
//                   
//
//----------------------------------------------------------------------

#include <string>
#include <iostream>
//#include <vector>
//#include <fstream>
//#include <getopt.h>
//#include <ctime>
//#include "Node.h"
#include "SemiflexibleGel.h"
//#include "Model.h"
//#include "Lbfgsb.h"
//#include "Lbfgs.h"
//#include "CGfast.h"
//#include "BrownianRod.h"
//#include <random/uniform.h>
//#include <random/normal.h>
////#include "Dirichlet.h"
//#include "PeriodicBox.h"
//#include "LeesEdwards.h"
//#include "GelOutput.h"
//#include "ViscousRegularizer.h"
//#include <process.h>
#include "SemiflexibleInput.h"
//#include "voom.h"
#include <fstream>
#include <cstdlib>
//#include <complex>
#include <map> // ?
#include <cmath> // maybe not needed
#include <sstream> //?
#include <iomanip> //?
#include <blitz/array.h> //?


using namespace std;
//using namespace tvmet;

//typedef std::map< std::string, std::string > ParamMap;

int main(int argc, char* argv[]){
	cout << "-----Start of program------" << endl;

	SemiflexibleInput inp(argv[1]);

	double numReal = -2.0;
	int numInt = -2.0;
	std::string numStr = "Blah blah";

	std::string testString = "kT";

	// Test accessors
	numReal = inp.getReal(testString);
	numInt = inp.getInt(testString);
	numStr = inp.getString(testString);

	cout << "\n" << "-----Accessor Testing------" << endl;
	cout << "This is numReal: " << numReal << endl;
	cout << "This is numInt: " << numInt << endl;
	cout << "This is numStr: " << numStr << endl;

	// Test mutators
	numReal = 3.14;
	numInt = 3.14;
	numStr = "3.14e0";	

	cout << "\n" << "-----Mutator Testing-------" << endl;
	inp.setReal(testString, numReal);
	numReal = inp.getReal(testString);
	cout << "This is numReal: " << numReal << endl;
	inp.setInt(testString, numInt);
	numInt = inp.getInt(testString);
	cout << "This is numInt: " << numInt << endl;
	inp.setStr(testString, numStr);
	numStr = inp.getString(testString);
	cout << "This is numStr: " << numStr << endl;
	// ParamMap propMap = inp.pm;


// for(ParamMap::const_iterator MapIterator = propMap.begin(); MapIterator != propMap.end(); ++MapIterator)
//   {
//    std::cout << "Key: \"" << MapIterator->first << "\" "
//    << "Value: " << MapIterator->second << endl;
//   }
// 	

	
	SemiflexibleGel<2> testInputGel (&inp);
	
	
	
	cout << "\n" << "-----More Testing----------" << endl;
	cout << "Bond Stiffness: " << inp.getReal("bond stiffness") << endl;
	cout << "Angle Stiffness: " << inp.getReal("angle stiffness") << endl;
	cout << "\n" << "-----End of program--------" << endl;
	
	
	return 0;
	}






// Functions

