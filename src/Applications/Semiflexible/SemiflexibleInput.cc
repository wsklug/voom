// -*- C++ -*-sint
//----------------------------------------------------------------------
//
//                          
//                
//                   
//
//----------------------------------------------------------------------


#include "SemiflexibleInput.h"


using namespace std;
using namespace tvmet;
using namespace voom;

//typedef SemiflexibleGel<2> Gel;
typedef std::map< std::string, std::string > ParamMap; 
typedef ParamMap::iterator PMIter;

#define USEOPENMP 0
#if defined(_OPENMP)
#include <omp.h>
#define USEOPENMP 1
#endif

    
// Verification -------------------------------------------
bool SemiflexibleInput::checkMap(std::string name) const
{
  if (_pm.find(name) == _pm.end()) {
		return false;
  }
  else {
  	return true;
  }
}




// Accessors ----------------------------------------------
double SemiflexibleInput::getReal(std::string name) const
{

  if (_pm.find(name) == _pm.end()) {
  	std::cerr << "Error: Value for " << name 
  	<< " not found in parameter map. \n"
  	<< "Check input file." << std::endl;
  	exit(1);
  }
  return atof((_pm.find(name)->second).c_str());
}

int SemiflexibleInput::getInt(std::string name) const
{
  if (_pm.find(name) == _pm.end()) {
    std::cerr << "Error: Value for " << name
	      << " not found in parameter map. \n"
	      << "Check input file." << std::endl;
    exit(1);
  }
  return atoi((_pm.find(name)->second).c_str());
}

std::string SemiflexibleInput::getString(std::string name) const
{
  if (_pm.find(name) == _pm.end()) {
    std::cerr << "Error: Value for " << name
	      << " not found in parameter map. \n"
	      << "Check input file." << std::endl;
    exit(1);
  }
  return _pm.find(name)->second;
}

bool SemiflexibleInput::getBool(std::string name) const
{
  if (_pm.find(name) == _pm.end()) {
    std::cerr << "Error: Value for " << name
	      << " not found in parameter map. \n"
	      << "Check input file." << std::endl;
    exit(1);
  }
  if (_pm.find(name)->second == "true")
    return true;
  else
    return false;	
}

// Mutators ----------------------------------------------- 
void SemiflexibleInput::setReal(std::string nameReal, double valueReal)
{
	std::ostringstream tmpVal;
	tmpVal << setprecision(16) << valueReal;
	_pm[nameReal] = tmpVal.str();
}	

void SemiflexibleInput::setInt(std::string nameInt, int valueInt)
{
	std::ostringstream tmpVal;
	tmpVal << valueInt;
	_pm[nameInt] = tmpVal.str();
}	

void SemiflexibleInput::setStr(std::string name, std::string value)
{
	_pm[name] = value;
}

void SemiflexibleInput::setBool(std::string name, bool value)
{
	_pm[name] = value;
}
// Previous Code ------------------------------------------
int getFilSegs_input(double avgCL, double nodesPerCL) {
  return (int)((avgCL+1)*nodesPerCL);
}


int getCurGelNum_input(std::string firstGelName) {  
  std::ifstream testStream;
  bool fileFound = true;
  int startPos = firstGelName.find("gelnum");
  int endPos = firstGelName.find(".gelsave");
  int i = atoi((firstGelName.substr(startPos+7,endPos-7-startPos)).data());
  
  // find first available gel name //
  while(fileFound) {
    testStream.open(firstGelName.c_str());
    if(testStream.fail()) {
      fileFound = false;
    }
    else{
      i++;
      startPos = firstGelName.find("gelnum");
      endPos = firstGelName.find(".gelsave");
      ostringstream tmpi;
      tmpi << i;
      firstGelName.replace(startPos+7,endPos-7-startPos,tmpi.str());
    }  
    testStream.close();
  }
  
  return i;
}

std::string getGelFileName_input(std::string gelLibDir, ParamMap & pmap, int gelNum=-1) {
  // get necessary parameters from map, construct file name //
  char * gelStoreFN = new char[256];
  double L = atof(pmap["L"].data());
  double l_c = L/(atof(pmap["L/l_c"].data()));
  double kcl = atof(pmap["crosslink stiffness"].data());
  double Wx = atof(pmap["Wx"].data());
  double Wy = atof(pmap["Wy"].data());
  double dL = atof(pmap["dL"].data());
  double Snem = atof(pmap["nematic order parameter"].data());
  std::string gelStoreName;
  
  if(gelNum == -1){
    sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-S=%f-gelnum=1.gelsave",L,l_c,dL,Wx,Wy,kcl,Snem);
    gelStoreName.assign(gelStoreFN);
    gelStoreName.insert(0,gelLibDir);
    gelNum = getCurGelNum_input(gelStoreName);
    gelStoreName.clear();
    delete gelStoreFN;
    gelStoreFN = new char[128];
  }
  
  sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-S=%f-gelnum=%d.gelsave",L,l_c,dL,Wx,Wy,kcl,Snem,gelNum);
  gelStoreName.assign(gelStoreFN);
  gelStoreName.insert(0,gelLibDir);
  return gelStoreName;
}

    ////////////////////////////////////////////////////////////////////
    // Constructor
    ////////////////////////////////////////////////////////////////////
SemiflexibleInput::SemiflexibleInput(std::string paramFileName)
{
    std::ifstream inFile(paramFileName.c_str());
    if(!(inFile.good())) {
        std::cerr << "Error: input file name does not exist!" << std::endl;
        exit(1);
    }
		std::string pfn(paramFileName);
    std::cout << "Input file name: " << pfn << std::endl;
    
    if(USEOPENMP) std::cout << "Open MP being used." << std::endl;
    else std::cout << "Open MP not being used." << std::endl;
    
    
    ////////////////////////////////////////////////////////////////////
    // Parameters
    ////////////////////////////////////////////////////////////////////
    
    double kT = -1.0; // temperature in pN-microns //
    double kC = -1.0; // bending modulus of actin in pN-microns^2 //
    double L_p = -1.0; // persistence length of actin in microns //
    double visc = -1.0; // viscosity //
    //  double r = -1.0; // radius of actin filaments in microns //
    double l_B = -1.0; // ratio of bending/stretching modulii //
    tvmet::Vector<double,2> syssize; // system size in microns //
    double kBond = -1.0; // spring constant for rods in pN/micron //
    double kAngle = -1.0; // spring constant for rod joints (angles) in pN-microns //
    double kcl = -1.0; // crosslink spring constant (negative for pinned node crosslinks) in pN/micron //
    double dt = -1.0; // time step in seconds //
    double L = -1.0; // length of filaments in microns //
    double dL = -1.0; // length of rods (filament sections) in microns //
    double L_over_lc = -1.0; // nominal value of L/l_c //
    double actuall_c = -1.0; // actual value of l_c //
    double lambda = -1.0; // l_c*(l_c/l_B)^z //
    double k_max = -1.0; // maximum stiffness allowed in nonlinear rods (entropic elasticity) //
    double entropic_k = -1.0;
    double entropic_Lp = -1.0;
    double filDens = -1.0; // filament density in microns^-2 //
    int nNodesPerFilament = -1; // # of nodes per filament //
    double nodesPerCL = 2.5; // avg # of nodes between crosslinks //
    double nematicOP = 0.0;
    double actualnemOP = -1.0;
    double nemDirectorAngle = 0.0;
    std::string nemPDFparam;
    bool prestress = false;
    std::string nemPDFtype = "Gaussian";
    
    int nGels2avg = 1;
    int curGelNum = -1;
    
    std::string gelDirectory = "./";
    std::string gelFileName;
    
    // JKP: This is weird. Should bondType be set as "Spring" in the input file 
    // instead of here? "bondType" is dependent on k_max, 
    // which is not listed in the input file. There is an F_max though.
    std::string bondType = "Spring";
    this->setStr("bondType",bondType);
    
    
    bool linearizedentropic = false;
    double linentfrac = 0.0;
    double linentmult = 1.0;
    
    std::string solverType = "LBFGSB";
    this->setStr("solverType","LBFGSB");
    bool cutOffEnds = false;
    this->setStr("Cut Off Ends","false");
    double viscReg = -1.0;
    
    bool relaxPrestress = false;
    this->setStr("relaxPrestress","false");
    bool retrieveGel = false;
    this->setStr("retrieveGel","false");
    bool checkConsist = false;
    this->setStr("checkConsist","false");
    bool shearXtest = false;
    this->setStr("ShearXTest","false");
    bool shearYtest = false;
    this->setStr("ShearYTest","false");
    bool expXYtest = false;
    this->setStr("expXYtest","false");
    bool expXtest = false;
    this->setStr("expXtest","false");
    bool expYtest = false;
    this->setStr("expYtest","false");
    
    bool writeStates = false;
    this->setStr("writeStates","false");
    bool zeroReturn = false;
    this->setStr("zeroReturn","false");
    bool adaptiveMeshing = false;
    this->setStr("adaptiveMeshing","false");
    double minLength = -1.0;
    
    std::string polydisp = "none";
    double actualL;

    
    ////////////////////////////////////////////////////////////////////
    // Motor Parameters
    ////////////////////////////////////////////////////////////////////
    
    double maxMotorForce = -1.0;
    double startMotorSep = -1.0;
    double annulusTol = -1.0;
    double motorDens = -1.0;
    
    ////////////////////////////////////////////////////////////////////
    // Nematic Test Parameters
    ////////////////////////////////////////////////////////////////////
    
    double shrStart = 0.0;
    _pm.insert(pair< std::string, std::string >("shrStart","false"));
    double shrEnd = -1.0;
    int nShrSteps = -1;
    double shrStep = -1.0;
    double stretchEnd = -1.0;
    int nStretchSteps = -1;
    double stretchStep = -1.0;
    double expandEnd = -1.0;
    int nExpandSteps = -1;
    double expandStep = -1.0;

    
    std::string curString;
    std::string parName;
    std::string parValStr;
    double parVal;
    std::string::iterator curStrIt;
    
    // Read input file
    
    while(!inFile.eof()) {
        char curLine[256];
        inFile.getline(curLine,256);
        curString.assign(curLine);
        curStrIt = curString.begin();
        if(*curStrIt != '#' && curString.length() > 2) {
            int pNameEnd = curString.find_first_of('\t');
            parName.assign(curString,0,pNameEnd);
            pNameEnd = curString.find_first_not_of('\t',pNameEnd);
            pNameEnd = curString.find_first_of('\t',pNameEnd);
            int pNameBegin = curString.find_first_not_of('\t',pNameEnd);
            pNameEnd = curString.find_first_of('\t',pNameBegin);
            parValStr.assign(curString,pNameBegin,pNameEnd-pNameBegin);
            if(parName.find("GelFile") != std::string::npos) {
	      gelFileName.assign(parValStr);
	      _pm.insert(pair< std::string, std::string>
			 ("storage file name",parValStr));
	    }
            else if(parName.find("GelDirectory") != std::string::npos){
	      gelDirectory.assign(parValStr);
	      _pm.insert(pair< std::string, std::string>
                         ("storage directory name",parValStr));
	    }
            else if(parName.find("kT")!=std::string::npos) {
	      kT = atof(parValStr.data());
	      _pm.insert(pair< std::string, std::string>("kT",parValStr));
            }
            else if(parName.find("L_p")!=std::string::npos) L_p = atof(parValStr.data());
            else if(parName.find("kC")!=string::npos) kC = atof(parValStr.data());
            else if(parName.find("visc. reg.")!=string::npos) {
	      viscReg = atof(parValStr.data());
            }
            else if(parName.find("visc")!=string::npos) {
	      visc = atof(parValStr.data());
	      _pm.insert(pair< std::string, std::string >("viscosity",parValStr));
	      
            }
            else if(parName.find("SolverType")!=string::npos){
	      solverType = parValStr.data();
	      this->setStr("solverType",solverType);
	    }
	    // else if(parName.find("r")!=string::npos) r = parVal;
            else if(parName.find("k_max")!=string::npos) {
              k_max = atof(parValStr.data());
              _pm.insert(pair< std::string, std::string >("k_max",parValStr));
            	bondType = "EntropicSpring";
              _pm["bondType"] = bondType;
            }
            else if(parName.find("Entropic_lin")!=string::npos) {
                if(atof(parValStr.data()) > 0.5) {
                    linearizedentropic = true;
                    _pm.insert(pair< std::string, std::string >("Entropic_lin_springs","1"));
                }
            }
            else if(parName.find("Entr_lin_frac")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("Entropic_lin_stiff_frac",parValStr.data()));
                linentfrac = atof(parValStr.data());
            }
            else if(parName.find("Entr_lin_mult")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("Entropic_lin_stiff_mult",parValStr.data()));
                linentmult = atof(parValStr.data());
            }
            else if(parName.find("Entropic_k")!=string::npos) {
                entropic_k = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("Entropic_k",parValStr));
            }
            else if(parName.find("Entropic_Lp")!=string::npos) {
                entropic_Lp = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("Entropic_Lp",parValStr));
            }
            else if(parName.find("NodesPerCL")!=string::npos) {
                nodesPerCL = atof(parValStr.data());
            }
            else if(parName.find("dt")!=string::npos) {
                dt = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("time step",parValStr));
            }
            else if(parName.find("L/l_C")!=string::npos) L_over_lc = atof(parValStr.data());
            else if(parName.find("Wx")!=string::npos) syssize[0] = atof(parValStr.data());
            else if(parName.find("Wy")!=string::npos) syssize[1] = atof(parValStr.data());
            else if(parName.find("kcl")!=string::npos) {
                kcl = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("crosslink stiffness",parValStr));
            }
            else if(parName.find("l_B")!=string::npos) {
                l_B = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("l_B",parValStr));
            }
            else if(parName.find("fit")!=string::npos) {
                //fitOrder = atoi(parValStr.data());
                _pm.insert(pair< std::string, std::string >("fit order",parValStr));
            }
            else if(parName.find("Shear_start")!=string::npos) shrStart = atof(parValStr.data());
            else if(parName.find("Shear_final")!=string::npos) shrEnd = atof(parValStr.data());
            else if(parName.find("Shear_size")!=string::npos) shrStep = atof(parValStr.data());
            else if(parName.find("Stretch_final")!=string::npos) stretchEnd = atof(parValStr.data());
            else if(parName.find("Stretch_steps")!=string::npos) nStretchSteps = atoi(parValStr.data());
            else if(parName.find("Expand_final")!=string::npos) expandEnd = atof(parValStr.data());
            else if(parName.find("Expand_steps")!=string::npos) nExpandSteps = atoi(parValStr.data());
            //       else if(parName.find("Affine_method")!=string::npos) affmethod.assign(parValStr);
            //       else if(parName.find("Affine_shear")!=string::npos) affshear = atof(parValStr.data());
            //       else if(parName.find("Affine_gels")!=string::npos) nGels2avg = atoi(parValStr.data());
            //       else if(parName.find("Affine_space")!=string::npos) affbinspace = atof(parValStr.data());
            //       else if(parName.find("Affine_tol")!=string::npos) affbintol = atof(parValStr.data());
            //       else if(parName.find("Affine_min")!=string::npos) affmin = atof(parValStr.data());
            //       else if(parName.find("Affine_max")!=string::npos) affmax = atof(parValStr.data());
            else if(parName.find("FilDens")!=string::npos) { 
            filDens = atof(parValStr.data());
            _pm.insert(pair< std::string, std::string >("filDens",parValStr.data()));
            }
            else if(parName.find("NemPDF_param")!=string::npos) {
                nemPDFparam = parValStr;
                _pm.insert(pair< std::string, std::string >("nematic PDF param",nemPDFparam));
            }
            else if(parName.find("NemPDF_type")!=string::npos) {
                nemPDFtype = parValStr;
                _pm.insert(pair< std::string, std::string >("orientational PDF",nemPDFtype));
            }
            else if(parName.find("NematicOP")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("nematic order parameter",parValStr));
                nematicOP = atof(parValStr.data());
            }
            else if(parName.find("NemDirAngle")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("nematic direction angle",parValStr));
                nemDirectorAngle = atof(parValStr.data());
            }
            else if(parName.find("MaxPrestress")!=string::npos) {
                if(atof(parValStr.data()) > 1.0e-6) {
                    _pm.insert(pair< std::string, std::string >("maximum prestress",parValStr));
                }
            }
            else if(parName.find("RelaxPrestress")!=string::npos) {
	      if(atof(parValStr.data()) >= .5){
		relaxPrestress = true;
		this->setStr("relaxPrestress","true");
	      }
            }
            else if(parName.find("Prestress")!=string::npos) {
                if(atof(parValStr.data()) < 0.5) _pm.insert(pair< std::string, std::string >("prestress","false"));
                else _pm.insert(pair< std::string, std::string >("prestress","true"));
            }
            else if(parName.find("Retrieve")!=string::npos) {
            	// JKP: why set retrieveGel = false if is initialized as false?
                if(atof(parValStr.data()) < 0.5) retrieveGel = false;
                else {retrieveGel = true;
		  //_pm["retrieveGel"] = "true";
		  this->setStr("retrieveGel","true");
                }
            }
            else if(parName.find("AdaptiveMesh")!=string::npos) {
                if(atof(parValStr.data()) >= .5) {
		  adaptiveMeshing = true;
		  //_pm["adaptiveMeshing"] = "true";
		  this->setStr("adaptiveMeshing","true");
		}
            }
            else if(parName.find("ShortSegRelief")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("nearby pair removal method",parValStr.data()));
            }
            else if(parName.find("TargetSegLength")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("target segment length", parValStr.data()));
            }
            else if(parName.find("MinSegLength")!=string::npos) {
                minLength = atof(parValStr.data()); // why is this here>?
                // minLength vs "minlength" is very confusing
                _pm.insert(pair< std::string, std::string >("Min Seg Length",
                parValStr.data()));
            }
            else if(parName.find("CutOffEnds")!=string::npos) {
	      if(atof(parValStr.data()) >= .5){
		cutOffEnds = true;
		_pm["Cut Off Ends"] = "true";
	      }
	    }
            else if(parName.find("ShearXTest")!=string::npos) {
	      if(atof(parValStr.data()) >= .5){
		shearXtest = true;
		_pm["ShearXTest"] = "true";
	      }
	    }
            else if(parName.find("ShearYTest")!=string::npos) {
                if(atof(parValStr.data()) >= .5) shearYtest = true;
            }
            else if(parName.find("ExpXYTest")!=string::npos) {
                if(atof(parValStr.data()) >= .5) expXYtest = true;
            }
            else if(parName.find("ExpXTest")!=string::npos) {
                if(atof(parValStr.data()) >= .5) expXtest = true;
            }
            else if(parName.find("ExpYTest")!=string::npos) {
                if(atof(parValStr.data()) >= .5) expYtest = true;
            }
            else if(parName.find("WriteStates")!=string::npos) {
	      if(atof(parValStr.data()) >= .5){
		writeStates = true;
		this->setStr("writeStates","true");
	      }
            }
            else if(parName.find("Zero_Return")!=string::npos) {
	      if(atof(parValStr.data()) >= .5){
		zeroReturn = true;
		this->setStr("zeroReturn","true");
	      }
	    }
            else if(parName.find("Polydispersity")!=string::npos) {
                polydisp = parValStr;
                _pm.insert(pair< std::string, std::string >("polydispersity",polydisp));
            }
            else if(parName.find("Long/short")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("longshortratio",parValStr));
            }
            else if(parName.find("Long fraction")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("longfraction",parValStr));
            }
            else if(parName.find("Long stiffness")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("longstiffness",parValStr));
            }
            else if(parName.find("Long bend_stiff")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("longbendstiffness",parValStr));
            }
            else if(parName.find("Long cutoff")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("longcutoff",parValStr));
            }
            else if(parName.find("Min. length")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("minlength",parValStr));
            }
            else if(parName.find("Check_consist")!=string::npos) {
                if(atof(parValStr.data()) >= .5) checkConsist = true;
            }
            else if(parName.find("Motor_dens")!=string::npos) {
                motorDens = atof(parValStr.data());
            }
            else if(parName.find("Max_Motor_force")!=string::npos) {
                maxMotorForce = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("motor force",parValStr));
            }
            else if(parName.find("Motor_sep_tol")!=string::npos) {
                annulusTol = atof(parValStr.data());
            }
            else if(parName.find("Motor_sep")!=string::npos) {
                startMotorSep = atof(parValStr.data());
            }
            else if(parName.find("L")!=string::npos) {
                L = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("L",parValStr));
            }
            // add polydispersity stuff //
        }

        curString.clear();
        parName.clear();
        parValStr.clear();
    }
    inFile.close();   

std::cout << "\n" << "Input file now closed." << "\n" << std::endl;
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
//                        Done reading input file
//
// Now we need to do a few calculations of "derived" parameters
// from those read from the input file.
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
 
// JKP: kC is calculated here? Why is it a commented out parameter in
// the input file?
//
// WSK: I think the bending modulus used to be a parameter that we
// could set in the input file, but then someone changed it to be
// derived from the persistence length and kT.
 if(kC < 0.0) {
   kC = kT*L_p;
 }	
 
 ///////////////////////////////////////////////////////////////////
 //
 // Retrieve gel geometry from a file
 //
 ///////////////////////////////////////////////////////////////////
 if(retrieveGel) {
      if(kC <= 0.0 || l_B <= 0.0) {
	std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
	exit(1);
      }
      
      // std::ostringstream tmpStr;
      // tmpStr << setprecision(16) << kC;
      if(!adaptiveMeshing) {

	// WSK: What is the difference between "bending modulus" and
	// "angle stiffness"?  If someone was using this as a clever
	// way to distinguish between adaptive and non-adaptive
	// meshing, that was foolish and should be avoided in favor of
	// some more explicit and clear strategy.
			this->setReal("bending modulus", kC); 
			}
			else {
			this->setReal("angle stiffness", kC); 
      }

      gelFileName.insert(0,gelDirectory);
      lambda = (L_over_lc/L)*pow(L_over_lc/(L*l_B),1.0/3.0);
      
      this->setReal("lambda", lambda);
      //input.setReal("lambda", lambda); 

      double mu = kC/sqrt(l_B);
      this->setReal("bond stiffness", mu);
      //input.setReal("bond stiffness", mu); 

    }

    ///////////////////////////////////////////////////////////////////
    //
    // Don't retrieve gel geometry from a file. 
    //
    ///////////////////////////////////////////////////////////////////
    else {

      //
      // WSK: This is wrong.  If the input file calls for anisotropy,
      // all of nematic order details should be listed there.
      //
      if(nematicOP >= 1.0e-6) {
			std::string nemPDFparam("Gaussian");
			_pm.insert(pair< std::string, std::string >("nematic PDFparam",nemPDFparam));
      }
    }
      //
      // WSK: This looks like the system size (of the box) is being
      // scaled by the filament length.  Is that the way it is
      // supposed to work?  That seems wacky. 
      //
      // At least the stored properties should be changed.
      //
      for(int dn=0;dn<2;dn++) {
			syssize[dn] *= L;
      }
    
      // get number of rods per filament from requirement of nodes
      // between crosslinks, then set all other parameters //
      if(!adaptiveMeshing) {
			nNodesPerFilament = 1 + getFilSegs_input(L_over_lc-1,nodesPerCL);
			dL = L/(nNodesPerFilament - 1);
			// Estimate kAngle = E*I/L = kT\xi_p/dL; kT=4.1pN-nm, \xi_p=10^4nm, 
			//dL~100nm
			//make kAngle 100 times larger to simulate the rotation diffusion
    
			if(kC > 0.0) {
	  	kAngle = kC/dL;
	  // std::ostringstream tmpStr;
	  // tmpStr << setprecision(16) << kAngle;
	  // _pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
	  this->setReal("angle stiffness", kAngle);
	}
	else {
	  std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
	  exit(1);
	}
	if(l_B > 0.0) {
	  kBond = kAngle/sqrt(l_B);
	  // std::ostringstream tmpStr;
	  // tmpStr << setprecision(16) << kBond;
	  // _pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr.str()));
	  this->setReal("bond stiffness", kBond);
	}
	else {
	  std::cerr << "Error: input file must have a value for the bending/stretching 			length l_B." << std::endl;
	}
      }
      //
      // Do adaptive meshing
      //
      else {
	dL = -1.0;
	// std::ostringstream tmpStr;
	// tmpStr << setprecision(16) << kC;
	// _pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
	
	//
	// With adaptive meshing the angle spring stiffness will not
	// be constant.  So in place of the angle stiffness, the
	// bending modulus is stored.  Why not just call it "bending
	// modulus"?
	//
	this->setReal("angle stiffness", kC);

	double mu = kC/sqrt(l_B);
	// std::ostringstream tmpStr2;
	// tmpStr2 << setprecision(16) << mu;
	// _pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr2.str()));

	// Same for bond stiffness; not constant when segment lengths
	// vary. Instead store stretching modulus.  Why not call it
	// "stretching modulus"?
	this->setReal("bond stiffness", mu);

      }

      lambda = (L_over_lc/L)*pow(L_over_lc/(L*l_B),1.0/3.0);

      this->setReal("lambda", lambda);
      
      this->setReal("Wx", syssize[0]);
      
      this->setReal("Wy", syssize[1]);
      
      this->setReal("L/l_c", L_over_lc);
      
      this->setReal("dL", dL);

      this->setReal("shrStart",shrStart);
      
      this->setReal("shrEnd",shrEnd);
      
      this->setReal("shrStep",shrStep);

      std::cout << "Constructing a gel with the following properties:" << std::endl;
      std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
      std::cout << "# nodes/filament: " << nNodesPerFilament << std::endl;
      std::cout << "Bond type: " << bondType << std::endl;
      //for(PMIter pmi=_pm.begin(); pmi!=_pm.end(); pmi++) {
      //std::cout << pmi->first << ": " << pmi->second << std::endl;
      //}

      std::string fName = getGelFileName_input(gelDirectory,_pm);
      
      curGelNum = getCurGelNum_input(fName);
      char curGelNumS [50];
      sprintf(curGelNumS, "%d",curGelNum);
      
      this->setStr("storage file name",fName);

      this->setInt("gel current number",curGelNum);
	
    
    // Write out parameter map
    std::cout << "\n" << "-----Parameter Map---------" << std::endl;
    
    for(ParamMap::const_iterator MapIterator = _pm.begin(); MapIterator != _pm.end(); ++MapIterator)
    {
        std::cout << "Key: \"" << MapIterator->first << "\" "
        << "Value: " << MapIterator->second << endl;
    }
}



