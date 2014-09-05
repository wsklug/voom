// -*- C++ -*-
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

typedef SemiflexibleGel<2> Gel;
typedef std::map< std::string, std::string > ParamMap; 
typedef ParamMap::iterator PMIter;

#define USEOPENMP 0
#if defined(_OPENMP)
#include <omp.h>
#define USEOPENMP 1
#endif

// Accessors ----------------------------------------------
double SemiflexibleInput::getReal(std::string name) const
{
  /* Cannot use map::operator[] to access value since it inserts a new element when no
  element matches the key */
  return atof((_pm.find(name)->second).c_str());
}

int SemiflexibleInput::getInt(std::string name) const
{
  return atoi((_pm.find(name)->second).c_str());
}

std::string SemiflexibleInput::getStr(std::string name) const
{
  return _pm.find(name)->second;
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

void SemiflexibleInput::setStr(std::string nameStr, std::string valueStr)
{
	_pm[nameStr] = valueStr;
}

// Previous Code ------------------------------------------
int getFilSegs(double avgCL, double nodesPerCL) {
  return (int)((avgCL+1)*nodesPerCL);
}


int getCurGelNum(std::string firstGelName) {  
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

std::string getGelFileName(std::string gelLibDir, ParamMap & pmap, int gelNum=-1) {
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
    gelNum = getCurGelNum(gelStoreName);
    gelStoreName.clear();
    delete gelStoreFN;
    gelStoreFN = new char[128];
  }
  
  sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-S=%f-gelnum=%d.gelsave",L,l_c,dL,Wx,Wy,kcl,Snem,gelNum);
  gelStoreName.assign(gelStoreFN);
  gelStoreName.insert(0,gelLibDir);
  return gelStoreName;
}

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
    double tmpratio = -1.0; // nominal value of L/l_c //
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
    std::string bondType = "Spring";
    
    bool linearizedentropic = false;
    double linentfrac = 0.0;
    double linentmult = 1.0;
    
    std::string solverType = "LBFGSB";
    
    bool cutOffEnds = false;
    
    double viscReg = -1.0;
    
    bool relaxPrestress = false;
    
    bool retrieveGel = false;
    
    bool checkConsist = false;
    
    bool shearXtest = false;
    bool shearYtest = false;
    bool expXYtest = false;
    bool expXtest = false;
    bool expYtest = false;
    
    bool writeStates = false;
    bool zeroReturn = false;
    
    bool adaptiveMeshing = false;
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
            if(parName.find("GelFile") != std::string::npos) gelFileName.assign(parValStr);
            else if(parName.find("GelDirectory") != std::string::npos) gelDirectory.assign(parValStr);
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
            else if(parName.find("SolverType")!=string::npos) solverType = parValStr.data();
            // else if(parName.find("r")!=string::npos) r = parVal;
            else if(parName.find("k_max")!=string::npos) {
                k_max = atof(parValStr.data());
                _pm.insert(pair< std::string, std::string >("k_max",parValStr));
                bondType = "EntropicSpring";
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
            else if(parName.find("L/l_C")!=string::npos) tmpratio = atof(parValStr.data());
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
            else if(parName.find("FilDens")!=string::npos) filDens = atof(parValStr.data());
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
                if(atof(parValStr.data()) >= .5) relaxPrestress = true;
            }
            else if(parName.find("Prestress")!=string::npos) {
                if(atof(parValStr.data()) < 0.5) _pm.insert(pair< std::string, std::string >("prestress","false"));
                else _pm.insert(pair< std::string, std::string >("prestress","true"));
            }
            else if(parName.find("Retrieve")!=string::npos) {
                if(atof(parValStr.data()) < 0.5) retrieveGel = false;
                else retrieveGel = true;
            }
            else if(parName.find("AdaptiveMesh")!=string::npos) {
                if(atof(parValStr.data()) >= .5) adaptiveMeshing = true;
            }
            else if(parName.find("ShortSegRelief")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("nearby pair removal method",parValStr.data()));
            }
            else if(parName.find("TargetSegLength")!=string::npos) {
                _pm.insert(pair< std::string, std::string >("target segment length", parValStr.data()));
            }
            else if(parName.find("MinSegLength")!=string::npos) {
                minLength = atof(parValStr.data());
            }
            else if(parName.find("CutOffEnds")!=string::npos) {
                if(atof(parValStr.data()) >= 0.5) cutOffEnds = true;
            }
            else if(parName.find("ShearXTest")!=string::npos) {
                if(atof(parValStr.data()) >= .5) shearXtest = true;
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
                if(atof(parValStr.data()) >= .5) writeStates = true;
            }
            else if(parName.find("Zero_Return")!=string::npos) {
                if(atof(parValStr.data()) >= .5) zeroReturn = true;
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

// -----------------------------------------
	  SemiflexibleGel<2> * gel;
    SemiflexibleGel<2>::DefNodeContainer nodes;
    PeriodicBox * box;
    GelOutput<2> output;
// -----------------------------------------

    // kC is calculated here? Why is it a commented out parameter in the input file?
    if(kC < 0.0) {
        kC = kT*L_p;
    }	
	
   if(retrieveGel) {
    if(kC > 0.0 && l_B > 0.0) {
      std::ostringstream tmpStr;
      tmpStr << setprecision(16) << kC;
      if(!adaptiveMeshing) _pm.insert(pair< std::string, std::string>("bending modulus",tmpStr.str()));
      else _pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
      // get nominal l_c and L from gel file name //
     //  int Lpos = gelFileName.find("L=") + 2;
//       int l_Cpos = gelFileName.find("l_C");
//       int dLpos = gelFileName.find("dL");
//       int Wxpos = gelFileName.find("Wx");
//       int kclpos = gelFileName.find("kcl");
//       int nempos = gelFileName.find("S");
//       int gelNumpos = gelFileName.find("gelnum");
//       int endPos = gelFileName.find(".gelsave");
//       std::string Lstr;
//       Lstr.assign(gelFileName,Lpos,l_Cpos-1-Lpos);
//       pm["L"] = Lstr;
//       L = atof(Lstr.data());
//       std::string l_Cstr;
//       l_Cstr.assign(gelFileName,l_Cpos+4,dLpos-5-l_Cpos);
//       tmpratio = L/atof(l_Cstr.data());
//       std::ostringstream tm;
//       tm << setprecision(16) << tmpratio;
//       pm["L/l_c"] = tm.str();
//       std::string dLstr;
//       dLstr.assign(gelFileName,dLpos+3,Wxpos-4-dLpos);
//       if(!adaptiveMeshing) {
//         pm["dL"] = dLstr;
//         dL = atof(dLstr.data());
//       }
//       else dL = -1.0;
//       std::string kclStr;
//       kclStr.assign(gelFileName,kclpos+4,nempos-5-kclpos);
//       kcl = atof(kclStr.data());
//       pm["crosslink stiffness"] = kclStr;
//       std::string nemstr;
//       nemstr.assign(gelFileName,nempos+2,gelNumpos-nempos-3);
//       nematicOP = atof(nemstr.data());
//       std::string GNstr;
//       GNstr.assign(gelFileName,gelNumpos+7,endPos-gelNumpos-7);
//       curGelNum = atoi(GNstr.data());
      
      gelFileName.insert(0,gelDirectory);
      lambda = (tmpratio/L)*pow(tmpratio/(L*l_B),1.0/3.0);
      std::ostringstream lambstr;
      lambstr << setprecision(16) << lambda;
      _pm["lambda"] = lambstr.str();
      
      double mu = kC/sqr(l_B);
      std::ostringstream mustream;
      mustream << setprecision(16) << mu;
      _pm["bond stiffness"] = mustream.str();
      
      
      if(!adaptiveMeshing) {
        gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutOffEnds,_pm,0.0);
	gel->compute(true,true,false);
	std::cout << "Sanity check: gel's energy at 0 shear = " << gel->energy() << std::endl;
      }
      else {
        gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutOffEnds,minLength,_pm);
      }
      box = gel->box();
      syssize = box->size();
      std::cout << "Retrieved and set up a gel with the following properties:" << std::endl;
      std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
      std::cout << "Bond type: " << bondType << std::endl;
      for(PMIter pmi=_pm.begin(); pmi!=_pm.end(); pmi++) {
	std::cout << pmi->first << ": " << pmi->second << std::endl;
      }
      std::ostringstream ss0;
      ss0 << setprecision(16) << syssize[0];
      std::ostringstream ss1;
      ss1 << setprecision(16) << syssize[1];
      _pm["Wx"] = ss0.str();
      _pm["Wy"] = ss1.str();
      
      if(relaxPrestress) {
	gel->removePrestress();
	std::cout << "Reset spring rest lengths/stiffnesses to relax prestress." << std::endl;
      }
    }
    else {
      std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
      exit(1);
    }
  }

  else {
    if(nematicOP >= 1.0e-6) {
      _pm.insert(pair< std::string, std::string >("nematic PDF param",nemPDFparam));
    }
    
    for(int dn=0;dn<2;dn++) {
      syssize[dn] *= L;
    }
    
    // get number of rods per filament from requirement of nodes between crosslinks, then set all other parameters //
    if(!adaptiveMeshing) {
      nNodesPerFilament = 1 + getFilSegs(tmpratio-1,nodesPerCL);
      dL = L/(nNodesPerFilament - 1);
      // Estimate kAngle = E*I/L = kT\xi_p/dL; kT=4.1pN-nm, \xi_p=10^4nm, dL~100nm
      //make kAngle 100 times larger to simulate the rotation diffusion
    
      if(kC > 0.0) {
        kAngle = kC/dL;
        std::ostringstream tmpStr;
        tmpStr << setprecision(16) << kAngle;
        _pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
      }
      else {
        std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
        exit(1);
      }
      if(l_B > 0.0) {
        kBond = kAngle/sqr(l_B);
        std::ostringstream tmpStr;
        tmpStr << setprecision(16) << kBond;
        _pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr.str()));
      }
      else {
        std::cerr << "Error: input file must have a value for the bending/stretching length l_B." << std::endl;
      }
    }
    else {
      dL = -1.0;
      std::ostringstream tmpStr;
      tmpStr << setprecision(16) << kC;
      _pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
      double mu = kC/sqr(l_B);
      std::ostringstream tmpStr2;
      tmpStr2 << setprecision(16) << mu;
      _pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr2.str()));
    }

    lambda = (tmpratio/L)*pow(tmpratio/(L*l_B),1.0/3.0);
    std::ostringstream lambstr;
    lambstr << setprecision(16) << lambda;
    _pm["lambda"] = lambstr.str();
    
    std::ostringstream tmpSx;
    tmpSx << setprecision(16) << syssize[0];
    _pm.insert(pair< std::string, std::string >("Wx",tmpSx.str()));
    std::ostringstream tmpSy;
    tmpSy << setprecision(16) << syssize[1];
    _pm.insert(pair< std::string, std::string >("Wy",tmpSy.str()));
    std::ostringstream tmpLl_c;
    tmpLl_c << setprecision(16) << tmpratio;
    _pm.insert(pair< std::string, std::string >("L/l_c",tmpLl_c.str()));
    std::ostringstream tmpdL;
    tmpdL << setprecision(16) << dL;
    _pm.insert(pair< std::string, std::string >("dL",tmpdL.str()));

    std::cout << "Constructing a gel with the following properties:" << std::endl;
    std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
    std::cout << "# nodes/filament: " << nNodesPerFilament << std::endl;
    std::cout << "Bond type: " << bondType << std::endl;
    for(PMIter pmi=_pm.begin(); pmi!=_pm.end(); pmi++) {
      std::cout << pmi->first << ": " << pmi->second << std::endl;
    }

    std::string fName = getGelFileName(gelDirectory,_pm);
    curGelNum = getCurGelNum(fName);

    if(adaptiveMeshing) _pm.insert(pair<std::string, std::string>("storage file name",fName));

    // make periodic box //
    box = new LeesEdwards(syssize[0],syssize[1],0.0);
    
    // create body //
    if(!adaptiveMeshing) {
      gel = new SemiflexibleGel<2>(nodes,box,filDens,nNodesPerFilament,dL,bondType,cutOffEnds,_pm);
      gel->compute(true,true,false);
      std::cout << "Sanity check: gel energy at zero shear = " << gel->energy() << std::endl;
    }
    else {
      gel = new SemiflexibleGel<2>(nodes,box,filDens,L,bondType,cutOffEnds,minLength,_pm);
    }
    //gel->addPinch(6.0,false,nodes,kBond,kAngle,visc,kT,dt,kcl);
    
    // write gel data to file //
    if(!adaptiveMeshing) {
      gel->storeGel(fName);
    }
}
    
    
    // Write out parameter map
    std::cout << "\n" << "-----Parameter Map---------" << std::endl;
    for(ParamMap::const_iterator MapIterator = _pm.begin(); MapIterator != _pm.end(); ++MapIterator)
    {
        std::cout << "Key: \"" << MapIterator->first << "\" "
        << "Value: " << MapIterator->second << endl;
    }


    
}
