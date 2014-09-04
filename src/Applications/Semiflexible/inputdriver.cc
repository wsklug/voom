// Test Input Driver for Semiflexible Class 

#include <string>
#include <iostream>
//#include <vector>
//#include <fstream>
//#include <getopt.h>
//#include <ctime>
//#include "Node.h"
//#include "SemiflexibleGel.h"
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

typedef std::map< std::string, std::string > ParamMap;

int main(){
    
    SemiflexibleInput inp("SFGel-nematictest12345.in");
	
	double numReal = -2.0;
	int numInt = -2.0;
	std::string numStr = "Blah blah";
	std::string okok = "kT";
	

	
	// Test accessors
	numReal = inp.getReal(okok);
	numInt = inp.getInt(okok);
	numStr = inp.getStr(okok);

	cout << "\n" << "-----Accessor Testing------" << endl;
	cout << "This is numReal: " << numReal << endl;
	cout << "This is numInt: " << numInt << endl;
	cout << "This is numStr: " << numStr << endl;
	
	// Test mutators
	numReal = 3.14;
	numInt = 3.14;
	numStr = "3.14e0";	

 	cout << "\n" << "-----Mutator Testing-------" << endl;
 	inp.setReal(okok, numReal);
 	numReal = inp.getReal(okok);
	cout << "This is numReal: " << numReal << endl;
	inp.setInt(okok, numInt);
	numInt = inp.getInt(okok);
 	cout << "This is numInt: " << numInt << endl;
 	inp.setStr(okok, numStr);
	numStr = inp.getStr(okok);
 	cout << "This is numStr: " << numStr << endl;
 	//ParamMap propMap = inp.pm;
	
	
// 	    for(ParamMap::const_iterator MapIterator = propMap.begin(); MapIterator != propMap.end(); ++MapIterator)
//     {
//         std::cout << "Key: \"" << MapIterator->first << "\" "
//         << "Value: " << MapIterator->second << endl;
//     }
// 	
    cout << "-----End of program-------" << endl;
    return 0;
}






// Functions

