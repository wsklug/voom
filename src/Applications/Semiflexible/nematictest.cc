#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include "Node.h"
#include "SemiflexibleGel.h"
#include "Model.h"
#include "Lbfgsb.h"
#include "Lbfgs.h"
#include "CGfast.h"
#include "BrownianRod.h"
#include <random/uniform.h>
#include <random/normal.h>
//#include "Dirichlet.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "GelOutput.h"
#include "ViscousRegularizer.h"
//#include <process.h>

#define USEOPENMP 0
#if defined(_OPENMP)
#include <omp.h>
#define USEOPENMP 1
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

typedef SemiflexibleGel<2> Gel;
typedef std::map< std::string, std::string > ParamMap;
typedef ParamMap::iterator PMIter;

struct NodePair {
  Gel::DefNode* node1;
  Gel::DefNode* node2;

  double initAngle;
};

struct DataPoint {
  double rL;

  int nDatPts;
  double mean;
  double stddev;
};

typedef std::vector< DataPoint > DataSet;
typedef std::vector< NodePair > NodePairList;
typedef std::vector< NodePairList > NodePairData;
typedef std::map< Gel::DefNode*, tvmet::Vector<double,2> > InitPositionMap;

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

// void setParamsfromFileName(std::string fileName, ParamMap & pmap) {
  
// }

void guessAffineShearX(Gel * gel, double shear) {
  Vector2D cent;
  cent = gel->box()->size();
  cent /= 2.0;
  int nFils = gel->filaments().size();
  for(int i=0; i<nFils; i++) {
    Gel::Filament * fil = gel->filament(i);
    int nNodesHere = fil->nodes.size();
    for(int j=0; j<nNodesHere; j++) {
      Vector2D diff;
      diff = fil->nodes[j]->position() - cent;
      Vector2D affdef;
      affdef[1] = 0.0;
      affdef[0] = shear*diff[1];
      affdef += fil->nodes[j]->position();
      fil->nodes[j]->setPoint(affdef);
    }
  }

}

void guessAffineShearY(Gel * gel, double shear) {
  Vector2D cent;
  cent = gel->box()->size();
  cent /= 2.0;
  int nFils = gel->filaments().size();
  for(int i=0; i<nFils; i++) {
    Gel::Filament * fil = gel->filament(i);
    int nNodesHere = fil->nodes.size();
    for(int j=0; j<nNodesHere; j++) {
      Vector2D diff;
      diff = fil->nodes[j]->position() - cent;
      Vector2D affdef;
      affdef[0] = 0.0;
      affdef[1] = -shear*diff[0];
      affdef += fil->nodes[j]->position();
      fil->nodes[j]->setPoint(affdef);
    }
  }

}

void guessExpXY(Gel * gel, double expand) {
  Vector2D cent;
  cent = gel->box()->size();
  cent /= 2.0;
  int nFils = gel->filaments().size();
  for(int i=0; i<nFils; i++) {
    Gel::Filament * fil = gel->filament(i);
    int nNodesHere = fil->nodes.size();
    for(int j=0; j<nNodesHere; j++) {
      Vector2D diff;
      diff = fil->nodes[j]->position() - cent;
      Vector2D affdef;
      affdef[1] = diff[1]*expand;
      affdef[0] = diff[0]*expand;
      affdef += fil->nodes[j]->position();
      fil->nodes[j]->setPoint(affdef);
    }
  }

}

void guessExpX(Gel * gel, double expand) {
  Vector2D cent;
  cent = gel->box()->size();
  cent /= 2.0;
  int nFils = gel->filaments().size();
  for(int i=0; i<nFils; i++) {
    Gel::Filament * fil = gel->filament(i);
    int nNodesHere = fil->nodes.size();
    for(int j=0; j<nNodesHere; j++) {
      Vector2D diff;
      diff = fil->nodes[j]->position() - cent;
      Vector2D affdef;
      affdef[1] = 0.0;
      affdef[0] = diff[0]*expand;
      affdef += fil->nodes[j]->position();
      fil->nodes[j]->setPoint(affdef);
    }
  }

}

void guessExpY(Gel * gel, double expand) {
  Vector2D cent;
  cent = gel->box()->size();
  cent /= 2.0;
  int nFils = gel->filaments().size();
  for(int i=0; i<nFils; i++) {
    Gel::Filament * fil = gel->filament(i);
    int nNodesHere = fil->nodes.size();
    for(int j=0; j<nNodesHere; j++) {
      Vector2D diff;
      diff = fil->nodes[j]->position() - cent;
      Vector2D affdef;
      affdef[0] = 0.0;
      affdef[1] = diff[1]*expand;
      affdef += fil->nodes[j]->position();
      fil->nodes[j]->setPoint(affdef);
    }
  }

}

int main(int argc, char* argv[]) {

  int verbose = 1;

  char* paramFileName = argv[1];

  std::ifstream inFile(paramFileName);
  if(!(inFile.good())) {
    std::cerr << "Error: input file name does not exist!" << std::endl;
    exit(1);
  }

  std::string pfn(paramFileName);
  int runnumStart = pfn.find_first_of("0123456789");
  int runnumEnd = pfn.find_first_not_of("0123456789",runnumStart);
  int runNum = atoi(pfn.substr(runnumStart,runnumEnd-runnumStart).data());

  char newDirName[128];
  sprintf(newDirName,"mkdir -p Run%d",runNum);
  system(newDirName);

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

  // motor parameters //
  double maxMotorForce = -1.0;
  double startMotorSep = -1.0;
  double annulusTol = -1.0;
  double motorDens = -1.0;

  // parameters for nematic test //
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
  
  ParamMap pm;
  std::string curString;
  std::string parName;
  std::string parValStr;
  double parVal;
  std::string::iterator curStrIt;
  
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
      if(parName.find("GelFile") != string::npos) gelFileName.assign(parValStr);
      else if(parName.find("GelDirectory") != string::npos) gelDirectory.assign(parValStr);
      else if(parName.find("kT")!=string::npos) {
	kT = atof(parValStr.data());
	pm.insert(pair< std::string, std::string>("kT",parValStr));
      }
      else if(parName.find("L_p")!=string::npos) L_p = atof(parValStr.data());
      else if(parName.find("kC")!=string::npos) kC = atof(parValStr.data());
      else if(parName.find("visc. reg.")!=string::npos) {
	viscReg = atof(parValStr.data());
      }
      else if(parName.find("visc")!=string::npos) {
	visc = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("viscosity",parValStr));
      }
      else if(parName.find("SolverType")!=string::npos) solverType = parValStr.data();
      // else if(parName.find("r")!=string::npos) r = parVal;
      else if(parName.find("k_max")!=string::npos) {
	k_max = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("k_max",parValStr));
	bondType = "EntropicSpring";
      }
      else if(parName.find("Entropic_lin")!=string::npos) {
	if(atof(parValStr.data()) > 0.5) {
	  linearizedentropic = true;
	  pm.insert(pair< std::string, std::string >("Entropic_lin_springs","1"));
	}
      }
      else if(parName.find("Entr_lin_frac")!=string::npos) {
	pm.insert(pair< std::string, std::string >("Entropic_lin_stiff_frac",parValStr.data()));
	linentfrac = atof(parValStr.data());
      }
      else if(parName.find("Entr_lin_mult")!=string::npos) {
	pm.insert(pair< std::string, std::string >("Entropic_lin_stiff_mult",parValStr.data()));
	linentmult = atof(parValStr.data());
      }
      else if(parName.find("Entropic_k")!=string::npos) {
	entropic_k = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("Entropic_k",parValStr));
      }
      else if(parName.find("Entropic_Lp")!=string::npos) {
	entropic_Lp = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("Entropic_Lp",parValStr));
      }
      else if(parName.find("NodesPerCL")!=string::npos) {
	nodesPerCL = atof(parValStr.data());
      }
      else if(parName.find("dt")!=string::npos) {
	dt = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("time step",parValStr));
      }
      else if(parName.find("L/l_C")!=string::npos) tmpratio = atof(parValStr.data());
      else if(parName.find("Wx")!=string::npos) syssize[0] = atof(parValStr.data());
      else if(parName.find("Wy")!=string::npos) syssize[1] = atof(parValStr.data());
      else if(parName.find("kcl")!=string::npos) {
	kcl = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("crosslink stiffness",parValStr));	
      }
      else if(parName.find("l_B")!=string::npos) {
	l_B = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("l_B",parValStr));
      }
      else if(parName.find("fit")!=string::npos) {
	//fitOrder = atoi(parValStr.data());
	pm.insert(pair< std::string, std::string >("fit order",parValStr));	
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
	pm.insert(pair< std::string, std::string >("nematic PDF param",nemPDFparam));
      }
      else if(parName.find("NemPDF_type")!=string::npos) {
	nemPDFtype = parValStr;
	pm.insert(pair< std::string, std::string >("orientational PDF",nemPDFtype));
      }
      else if(parName.find("NematicOP")!=string::npos) {
	pm.insert(pair< std::string, std::string >("nematic order parameter",parValStr));
	nematicOP = atof(parValStr.data());
      }
      else if(parName.find("NemDirAngle")!=string::npos) {
        pm.insert(pair< std::string, std::string >("nematic direction angle",parValStr));
	nemDirectorAngle = atof(parValStr.data());
      }
      else if(parName.find("MaxPrestress")!=string::npos) {
	if(atof(parValStr.data()) > 1.0e-6) {
	  pm.insert(pair< std::string, std::string >("maximum prestress",parValStr));
        }
      }
      else if(parName.find("RelaxPrestress")!=string::npos) {
	if(atof(parValStr.data()) >= .5) relaxPrestress = true;
      }
      else if(parName.find("Prestress")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) pm.insert(pair< std::string, std::string >("prestress","false"));
	else pm.insert(pair< std::string, std::string >("prestress","true"));
      }
      else if(parName.find("Retrieve")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) retrieveGel = false;
	else retrieveGel = true;
      }
      else if(parName.find("AdaptiveMesh")!=string::npos) {
	if(atof(parValStr.data()) >= .5) adaptiveMeshing = true;
      }
      else if(parName.find("ShortSegRelief")!=string::npos) {
	pm.insert(pair< std::string, std::string >("nearby pair removal method",parValStr.data()));
      }
      else if(parName.find("TargetSegLength")!=string::npos) {
	pm.insert(pair< std::string, std::string >("target segment length", parValStr.data()));
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
	pm.insert(pair< std::string, std::string >("polydispersity",polydisp));
      }
      else if(parName.find("Long/short")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longshortratio",parValStr));
      }
      else if(parName.find("Long fraction")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longfraction",parValStr));
      }
      else if(parName.find("Long stiffness")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longstiffness",parValStr));
      }
      else if(parName.find("Long bend_stiff")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longbendstiffness",parValStr));
      }
      else if(parName.find("Long cutoff")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longcutoff",parValStr));
      }
      else if(parName.find("Min. length")!=string::npos) {
	pm.insert(pair< std::string, std::string >("minlength",parValStr));
      }
      else if(parName.find("Check_consist")!=string::npos) {
	if(atof(parValStr.data()) >= .5) checkConsist = true;
      }
      else if(parName.find("Motor_dens")!=string::npos) {
	motorDens = atof(parValStr.data());
      }
      else if(parName.find("Max_Motor_force")!=string::npos) {
	maxMotorForce = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("motor force",parValStr));
      }
      else if(parName.find("Motor_sep_tol")!=string::npos) {
	annulusTol = atof(parValStr.data());
      }
      else if(parName.find("Motor_sep")!=string::npos) {
	startMotorSep = atof(parValStr.data());
      }
      else if(parName.find("L")!=string::npos) {
	L = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("L",parValStr));
      }
      // add polydispersity stuff //
    }
    curString.clear();
    parName.clear();
    parValStr.clear();
  }
  inFile.close();

  SemiflexibleGel<2> * gel;
  SemiflexibleGel<2>::DefNodeContainer nodes;
  PeriodicBox * box;
  GelOutput<2> output;
  
  if(kC < 0.0) {
    kC = kT*L_p;
  }

  if(retrieveGel) {
    if(kC > 0.0 && l_B > 0.0) {
      std::ostringstream tmpStr;
      tmpStr << setprecision(16) << kC;
      if(!adaptiveMeshing) pm.insert(pair< std::string, std::string>("bending modulus",tmpStr.str()));
      else pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
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
      pm["lambda"] = lambstr.str();
      
      double mu = kC/sqr(l_B);
      std::ostringstream mustream;
      mustream << setprecision(16) << mu;
      pm["bond stiffness"] = mustream.str();
      
      
      if(!adaptiveMeshing) {
        gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutOffEnds,pm,0.0);
	gel->compute(true,true,false);
	std::cout << "Sanity check: gel's energy at 0 shear = " << gel->energy() << std::endl;
      }
      else {
        gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutOffEnds,minLength,pm);
      }
      box = gel->box();
      syssize = box->size();
      std::cout << "Retrieved and set up a gel with the following properties:" << std::endl;
      std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
      std::cout << "Bond type: " << bondType << std::endl;
      for(PMIter pmi=pm.begin(); pmi!=pm.end(); pmi++) {
	std::cout << pmi->first << ": " << pmi->second << std::endl;
      }
      std::ostringstream ss0;
      ss0 << setprecision(16) << syssize[0];
      std::ostringstream ss1;
      ss1 << setprecision(16) << syssize[1];
      pm["Wx"] = ss0.str();
      pm["Wy"] = ss1.str();
      
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

  //
  // Don't retrieve gel, create it.
  //
  else {
    if(nematicOP >= 1.0e-6) {
      pm.insert(pair< std::string, std::string >("nematic PDF param",nemPDFparam));
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
        pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
      }
      else {
        std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
        exit(1);
      }
      if(l_B > 0.0) {
        kBond = kAngle/sqr(l_B);
        std::ostringstream tmpStr;
        tmpStr << setprecision(16) << kBond;
        pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr.str()));
      }
      else {
        std::cerr << "Error: input file must have a value for the bending/stretching length l_B." << std::endl;
      }
    }
    else {
      dL = -1.0;
      std::ostringstream tmpStr;
      tmpStr << setprecision(16) << kC;
      pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
      double mu = kC/sqr(l_B);
      std::ostringstream tmpStr2;
      tmpStr2 << setprecision(16) << mu;
      pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr2.str()));
    }

    lambda = (tmpratio/L)*pow(tmpratio/(L*l_B),1.0/3.0);
    std::ostringstream lambstr;
    lambstr << setprecision(16) << lambda;
    pm["lambda"] = lambstr.str();
    
    std::ostringstream tmpSx;
    tmpSx << setprecision(16) << syssize[0];
    pm.insert(pair< std::string, std::string >("Wx",tmpSx.str()));
    std::ostringstream tmpSy;
    tmpSy << setprecision(16) << syssize[1];
    pm.insert(pair< std::string, std::string >("Wy",tmpSy.str()));
    std::ostringstream tmpLl_c;
    tmpLl_c << setprecision(16) << tmpratio;
    pm.insert(pair< std::string, std::string >("L/l_c",tmpLl_c.str()));
    std::ostringstream tmpdL;
    tmpdL << setprecision(16) << dL;
    pm.insert(pair< std::string, std::string >("dL",tmpdL.str()));

    std::cout << "Constructing a gel with the following properties:" << std::endl;
    std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
    std::cout << "# nodes/filament: " << nNodesPerFilament << std::endl;
    std::cout << "Bond type: " << bondType << std::endl;
    for(PMIter pmi=pm.begin(); pmi!=pm.end(); pmi++) {
      std::cout << pmi->first << ": " << pmi->second << std::endl;
    }

    std::string fName = getGelFileName(gelDirectory,pm);
    curGelNum = getCurGelNum(fName);

    if(adaptiveMeshing) pm.insert(pair<std::string, std::string>("storage file name",fName));

    // make periodic box //
    box = new LeesEdwards(syssize[0],syssize[1],0.0);
    
    // create body //
    if(!adaptiveMeshing) {
      gel = new SemiflexibleGel<2>(nodes,box,filDens,nNodesPerFilament,dL,bondType,cutOffEnds,pm);
      gel->compute(true,true,false);
      std::cout << "Sanity check: gel energy at zero shear = " << gel->energy() << std::endl;
    }
    else {
      gel = new SemiflexibleGel<2>(nodes,box,filDens,L,bondType,cutOffEnds,minLength,pm);
    }
    //gel->addPinch(6.0,false,nodes,kBond,kAngle,visc,kT,dt,kcl);
    
    // write gel data to file //
    if(!adaptiveMeshing) {
      gel->storeGel(fName);
    }

    // now retrieve gel and confirm that nodes are at same positions //
//     SemiflexibleGel<2> * gel2;
//     SemiflexibleGel<2>::DefNodeContainer nodes2;
//     PeriodicBox * box2;
//     std::ostringstream tmpStr;
//     tmpStr << setprecision(16) << kC;
//     pm.insert(pair< std::string, std::string>("bending modulus",tmpStr.str()));
//     std::ostringstream kclStr;
//     kclStr << setprecision(16) << kcl;
//     pm["crosslink stiffness"] = kclStr.str();
//     gel2 = new SemiflexibleGel<2>(fName,nodes2,bondType,cutOffEnds,pm);
//     int nF = gel2->filaments().size();
//     for(int fi=0; fi<nF; fi++) {
//       SemiflexibleGel<2>::Filament * f1p = gel->filament(fi);
//       SemiflexibleGel<2>::Filament * f2p = gel2->filament(fi);
//       int nA2,nB2,nN2;
//       nA2 = f2p->angles.size();
//       nB2 = f2p->bonds.size();
//       nN2 = f2p->nodes.size();
//       int nA,nB,nN;
//       nA = f1p->angles.size();
//       nB = f1p->bonds.size();
//       nN = f1p->nodes.size();
//       if(nA!=nA2 || nB!=nB2 || nN!=nN2) {
// 	std::cout << "Error in gel storage/retrieval: filament " << fi << " was not stored or retrieved correctly." << std::endl;
//       }
//       else {
// 	SemiflexibleGel<2>::AngleContainer angs1;
// 	SemiflexibleGel<2>::AngleContainer angs2;
// 	angs1 = f1p->angles;
// 	angs2 = f2p->angles;
// 	for(int ai=0; ai<nA; ai++) {
// 	  if(angs1[ai]->stiffness() != angs2[ai]->stiffness()) {
// 	    std::cout << "Error in gel storage/retrieval: angle spring " << ai << " on filament " << fi << "." << std::endl;
// 	    std::cout << "Original stiffness = " << angs1[ai]->stiffness() << "; retrieved stiffness = " << angs2[ai]->stiffness() << "." << std::endl;
// 	  }
// 	}
// 	for(int bi=0; bi<nB; bi++) {
// 	  if(f1p->bonds[bi]->stiffness() != f2p->bonds[bi]->stiffness()) {
// 	    std::cout << "Error in gel storage/retrieval: spring " << bi << " on filament " << fi << "." << std::endl;
// 	    std::cout << "Original stiffness = " << f1p->bonds[bi]->stiffness() << "; retrieved stiffness = " << f2p->bonds[bi]->stiffness() << "." << std::endl;
// 	    std::cout << "Stiffness disparity = " << abs(f1p->bonds[bi]->stiffness()-f2p->bonds[bi]->stiffness()) << std::endl;
// 	  }
// 	}
// 	for(int ni=0; ni<nN; ni++) {
// 	  if(f1p->nodes[ni]->point()[0] != f2p->nodes[ni]->point()[0] || f1p->nodes[ni]->point()[1] != f2p->nodes[ni]->point()[1]) {
// 	    std::cout << "Error in gel storage/retrieval: node " << ni << " on filament " << fi << " has inconsistent displacement value." << std::endl;
// 	  }
// 	  if(f1p->nodes[ni]->position()[0] != f2p->nodes[ni]->position()[0] || f1p->nodes[ni]->position()[1] != f2p->nodes[ni]->position()[1]) {
// 	    std::cout << "Error in gel storage/retrieval: node " << ni << " on filament " << fi << " has inconsistent position value." << std::endl;
// 	  }
// 	}
//       }
//     }
    fName.clear();

    if(relaxPrestress) {
      gel->removePrestress();
      std::cout << "Reset spring rest lengths/stiffnesses to relax prestress." << std::endl;
    }

    char clinkdistfile[128];
    sprintf(clinkdistfile,"Run%d/CrosslinkSepDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,L,L/tmpratio,dL,nematicOP,curGelNum);
    fName.assign(clinkdistfile);
    output.printCrosslinkData(gel,fName);

    fName.clear();
    char fillendistfile[128];
    sprintf(fillendistfile,"Run%d/FilLengthDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,L,L/tmpratio,dL,nematicOP,curGelNum);
    fName.assign(fillendistfile);
    output.printFilLengthData(gel,fName);

    fName.clear();
    char nematicdistfile[128];
    sprintf(nematicdistfile,"Run%d/NematicDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,L,L/tmpratio,dL,nematicOP,curGelNum);
    fName.assign(nematicdistfile);
    output.printNematicData(gel,fName);
 
    fName.clear();
  }

  // check to make sure filaments aren't multiply crosslinked //
  // gel->checkCrosslinks();

  if(motorDens > 0.0 && adaptiveMeshing) {
    gel->addPinches(motorDens,startMotorSep,annulusTol,maxMotorForce,false);
  }

  if(linearizedentropic) {
    // print out list of stiffened segments //
    char stiffFileName[256];
    sprintf(stiffFileName,"Run%d/stiffsegs-%d.dat",runNum,runNum);
    std::string stiffFN(stiffFileName);
    gel->printStiffenedSegments(stiffFN,0.0);
  }

  if(shearXtest || shearYtest || expXYtest || expXtest || expYtest) {
    actuall_c = gel->getMeanCLsep();
    actualnemOP = gel->getNematicOP();
    double actuallamb = actuall_c*pow(actuall_c/l_B,1.0/3.0);

    char nemFileName[256];
    
    sprintf(nemFileName,"Run%d/nematicgel-L=%f-l_C=%f-dL=%f-S=%f-Wx=%f-Wy=%f-gelnum=%d",runNum,L,L/tmpratio,dL,nematicOP,syssize[0],syssize[1],curGelNum);
    
    //first, let system relax with no shear, stretch, etc., and store resultant positions in map.  then do various tests, returning after each to the original relaxed state//
    Model::NodeContainer modelnodes;
    modelnodes.reserve( nodes.size() );
    for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      modelnodes.push_back( *n );
    }
    Model model(modelnodes);
    model.pushBackBody( gel );
    
    int m=7;
    double factr=1.0e+1;
    double pgtol=1.0e-6;
    double gtol = 1.0e-5;
    int iprint = 1;
    int maxiter = 1500000;
    Solver * solver;
    if(solverType.find("LBFGSB")!=string::npos) {
      ifstream lbfgsbinp("lbfgsb.inp");
      lbfgsbinp >> iprint >> factr >> pgtol >> m ;
      if(verbose) {
	std::cout << "Input iprint: " << iprint << std::endl
		  << "Input factr: " << factr << std::endl
		  << "Input pgtol: " << pgtol << std::endl
		  << "Input m: " << m << std::endl;
      }

    }
    else {
      ifstream lbfgsinp("lbfgs.inp");
      lbfgsinp >> iprint >> gtol >> m ;
      if(verbose) {
	std::cout << "Input iprint: " << iprint << std::endl
		  << "Input gtol: " << gtol << std::endl
		  << "Input m: " << m << std::endl;
      }
 
    }

    ViscousRegularizer* vReg;
    double vRegTol = nodes.size()*viscReg*(1.0e-12);
    if(viscReg > 0.0) {
      gel->setViscReg(viscReg);
      vReg = gel->viscReg();
      ifstream viscreginp("viscreg.inp");
      viscreginp >> vRegTol >> maxiter;

      std::cout << "Viscous regularizer parameters: " << std::endl
		<< "Input visc. reg. energy tol: " << vRegTol << std::endl
		<< "Input solver iters between resets: " << maxiter << std::endl;
    }

    if(solverType.find("LBFGSB")!=string::npos) {
      solver = new Lbfgsb(model.dof(),m,factr,pgtol,iprint,maxiter);
    }

    else {
      solver = new Lbfgs(m,gtol,iprint,maxiter);
    }

    //    Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxiter );
//     box->setShear(0.0);
//     std::cout << "Relaxing gel with no shear/stretching." << std::endl;
//     solver->solve(&model);

    InitPositionMap ipmap;
    for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      ipmap[*n] = (*n)->point();
    }



    if(zeroReturn) {
      box->setShear(.0025);
      double oldenergy = gel->energy();
      std::cout << "Current energy with no shear = " << oldenergy << "; shearing gel slightly and relaxing, then setting shear to 0 and relaxing again." << std::endl;	     
      solver->solve(&model);
      box->setShear(0.0);
      solver->solve(&model);
      std::cout << "New energy with no shear = " << gel->energy() << "." << std::endl;
      if(gel->energy() <= oldenergy) {    
        for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
          ipmap[*n] = (*n)->point();
        }
      }
    }

    std::string fname2print;

    std::string nematicFileName(nemFileName);

//     char solverDataFN[256];
//     sprintf(solverDataFN,"Run%d/solverData",runNum);
//     std::string solverDataFileName(solverDataFN);
    std::ostringstream paramstring;
    paramstring << "Parameters: "
		<< "L/l_C = " << L/actuall_c
		<< ", L = " << L
		<< ", S = " << actualnemOP
		<< ", l_B = " << l_B
		<< ", Wx = " << syssize[0]
		<< ", Wy = " << syssize[1];
    

    // rotate system and make sure energy is the same //
//     if(abs(syssize[0]-syssize[1]) < 1.0e-6) {
//       double olden = gel->energy();
//       for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
// 	tvmet::Vector<double,2> oldpos = (*n)->point();
// 	tvmet::Vector<double,2> newpos;
// 	newpos[0] = -oldpos[1];
// 	newpos[1] = oldpos[0];
// 	(*n)->setPoint(newpos);
//       }
//       double newen = gel->energy();
//       solver.solve(&model);
//       std::cout << "Gel energy before rotation = " << olden << "; gel energy after rotation = " << newen << "; gel energy after new minimization = " << gel->energy() << "." << std::endl;
//     }


    if(checkConsist) {
      double perturbScale = minLength/100.0;
      //double perturbScale = 0.0;
      ranlib::Normal<double> perturbRNG(0.0,perturbScale);
      perturbRNG.seed((unsigned int)time(0));
      for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
	tvmet::Vector<double,2> oldPoint = (*n)->point();
	for(int ido=0; ido<2; ido++) oldPoint[ido] += perturbRNG.random();
	(*n)->setPoint(oldPoint);
      }
      model.checkConsistency(true,false);
    }

    if(shearXtest) {

      std::cout << std::endl << "Shearing gel in X direction." << std::endl;
      std::string shearXFileName = nematicFileName + "-ShearX.dat";
      std::string shearXAffineFileName = nematicFileName + "-affine-ShearX.dat";
      for(InitPositionMap::iterator ipit=ipmap.begin(); ipit!=ipmap.end(); ipit++) {
	(ipit->first)->setPoint(ipit->second);
      }
      output.printParamHeader(shearXFileName,paramstring.str());
      output.printParamHeader(shearXAffineFileName,paramstring.str());
      std::string fieldString = "ShearX\tE_tot\tE_gel\tE_bend\tE_strc\tE_cl\tE_mot\tE_pnch\tE_vreg";
      output.printFieldLabels(shearXFileName,fieldString);
      output.printFieldLabels(shearXAffineFileName,fieldString);
      nShrSteps = (int)(abs(shrEnd-shrStart)/shrStep);
      for(int sx=0; sx<=nShrSteps; sx++) {
	//gel->turnOffPinches();
	std::cout << "Shear X = " << shrStart + shrStep*sx << std::endl;
	box->setShearX(shrStart + shrStep*sx);
	guessAffineShearX(gel,shrStart + shrStep*sx);
	gel->compute(true,true,false);
	output.printEnergies(gel,shearXAffineFileName,shrStep*sx);
	if(viscReg > 0.0) {
	  vReg->step();
	  solver->solve(&model);
	  vReg->compute(true,true,false);
	  double viscEng = vReg->energy();
	  std::cout << "Viscous regularizer energy = " << viscEng << std::endl;
	  while(viscEng > vRegTol) {
	    vReg->step();
	    solver->solve(&model);
	    vReg->compute(true,true,false);
	    viscEng = vReg->energy();
	    std::cout << "Viscous regularizer energy = " << viscEng << std::endl;
	  }
	  vReg->step();
	}
	else {
	  solver->solve(&model);
	}
	gel->compute(true,true,false);
	output.printEnergies(gel,shearXFileName,shrStart+(shrStep*sx));
	char affXname[128];
	sprintf(affXname,"Run%d/afftest-shearX=%f.dat",runNum,shrStart+(shrStep*sx));
	std::string affXFN(affXname);
	gel->storeGel(affXFN);
	//if(motorDens > 0.0) {
	//  double motForceStep = maxMotorForce/5.0;
	//  for(int mi=1; mi<=5; mi++) {
	//    double motForce = motForceStep*mi;
	//    std::cout << "Motor force = " << motForce << "." << std::endl;
	//    gel->turnOnPinches(motForce);
	//    solver->solve(&model);
	//  }
	//  gel->compute(true,true,false);
	//  //gel->printBigBends();
	//  output.printEnergies(gel,shearXFileName,shrStart+(shrStep*sx));

	  //char affXname[128];
	  //sprintf(affXname,"Run%d/afftest-motors-shearX=%f.dat",runNum,shrStart+(shrStep*sx));
	  //std::string affXFN(affXname);
	  //gel->storeGel(affXFN);
	//}

	//output.printSolverData(gel,solver,shearXSolverFN);
	if(writeStates) {
	  fname2print.clear();
	  char shearXname[128];
	  sprintf(shearXname,"Run%d/shearX=%f",runNum,shrStart+(shrStep*sx));
	  fname2print.assign(shearXname);
	  //gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}
        
      }
    }

    if(shearYtest) {
      std::cout << std::endl << "Shearing gel in Y direction." << std::endl;
      for(InitPositionMap::iterator ipit=ipmap.begin(); ipit!=ipmap.end(); ipit++) {
	(ipit->first)->setPoint(ipit->second);
      }
      
      box->setShear(0.0);

      std::string shearYFileName = nematicFileName + "-ShearY.dat";
      std::string shearYAffineFileName = nematicFileName + "-affine-ShearY.dat";
      output.printParamHeader(shearYFileName,paramstring.str());
      output.printParamHeader(shearYAffineFileName,paramstring.str());	
      std::string fieldString = "ShearY\tE_tot\tE_gel\tE_bend\tE_strc\tE_cl\tE_mot\tE_pnch";
      //std::string fieldString = "ShearY\tE_tot\tE_cl\tE_bend\tE_strc";
      output.printFieldLabels(shearYFileName,fieldString);
      output.printFieldLabels(shearYAffineFileName,fieldString);
      nShrSteps = (int)(abs(shrEnd)/shrStep);
      for(int sy=0; sy<=nShrSteps; sy++) {
	std::cout << "Shear Y = " << shrStep*sy << std::endl;
	box->setShearY(shrStep*sy);
	guessAffineShearY(gel,shrStep*sy);
	gel->compute(true,true,false);
	output.printEnergies(gel,shearYAffineFileName,shrStep*sy);
	solver->solve(&model);
	gel->compute(true,true,false);
	output.printEnergies(gel,shearYFileName,shrStep*sy);
	//output.printSolverData(gel,solver,shearYSolverFN);
	if(writeStates) {
	  fname2print.clear();
	  char shearYname[128];
	  sprintf(shearYname,"Run%d/shearY=%f",runNum,shrStep*sy);
	  fname2print.assign(shearYname);
// 	  gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}

        char affYname[128];
	sprintf(affYname,"Run%d/afftest-shearY=%f.dat",runNum,shrStep*sy);
        std::string affYFN(affYname);
        gel->storeGel(affYFN);
      }
    }
    
    if(expXYtest) {
      std::cout << "Expanding gel uniformly." << std::endl;
      for(InitPositionMap::iterator ipit=ipmap.begin(); ipit!=ipmap.end(); ipit++) {
	(ipit->first)->setPoint(ipit->second);
      }

      box->setShear(0.0);

      std::string expXYFileName = nematicFileName + "-ExpXY.dat";
      std::string expXYAffineFileName = nematicFileName + "-affine-ExpXY.dat";
      output.printParamHeader(expXYFileName,paramstring.str());
      output.printParamHeader(expXYAffineFileName,paramstring.str());
      //std::string fieldString = "ExpXY\tE_tot\tE_cl\tE_bend\tE_strc";
      std::string fieldString = "ExpXY\tE_tot\tE_gel\tE_bend\tE_strc\tE_cl\tE_mot\tE_pnch";
      output.printFieldLabels(expXYFileName,fieldString);
      output.printFieldLabels(expXYAffineFileName,fieldString);
      
      expandStep = (expandEnd-1.0)/nExpandSteps;
      for(int ne=0; ne<=nExpandSteps; ne++) {
	double strtchFactor = 1.0+(expandStep*ne);
	double strtchPct = 100.0*(strtchFactor-1.0);
	std::cout << "Isotropic expansion pct. = " << strtchPct << std::endl;
	box->setSize(syssize);
	box->stretchXY(strtchFactor);
	guessExpXY(gel,expandStep*ne);
	gel->compute(true,true,false);
	output.printEnergies(gel,expXYAffineFileName,strtchFactor-1.0);
	solver->solve(&model);
	gel->compute(true,true,false);
	output.printEnergies(gel,expXYFileName,strtchFactor-1.0);
	//output.printSolverData(gel,solver,expXYSolverFN);
	if(writeStates) {
	  fname2print.clear();
	  char expXYname[128];
	  sprintf(expXYname,"Run%d/expXY=%f",runNum,strtchFactor);
	  fname2print.assign(expXYname);
	  //gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}

	char affExpXYname[128];
	sprintf(affExpXYname,"Run%d/afftest-expXY=%f.dat",runNum,strtchFactor);
        std::string affExpXYFN(affExpXYname);
        gel->storeGel(affExpXYFN);
      }
    }
    
    if(expXtest) {
      std::cout << "Stretching gel in X direction." << std::endl;
      for(InitPositionMap::iterator ipit=ipmap.begin(); ipit!=ipmap.end(); ipit++) {
	(ipit->first)->setPoint(ipit->second);
      }
      
      box->setShear(0.0);
      box->setSize(syssize);
      
      std::string expXFileName = nematicFileName + "-ExpX.dat";
      std::string expXAffineFileName = nematicFileName + "-affine-ExpX.dat";
      output.printParamHeader(expXFileName,paramstring.str());
      output.printParamHeader(expXAffineFileName,paramstring.str());
      //std::string fieldString = "StrX\tE_tot\tE_cl\tE_bend\tE_strc";
      std::string fieldString = "StrX\tE_tot\tE_gel\tE_bend\tE_strc\tE_cl\tE_mot\tE_pnch";
      output.printFieldLabels(expXFileName,fieldString);
      output.printFieldLabels(expXAffineFileName,fieldString);
      
      stretchStep = (stretchEnd-1.0)/nStretchSteps;
      for(int nsx=0; nsx<=nStretchSteps; nsx++) {
	double strtchFactor = 1.0+(stretchStep*nsx);
	double strtchPct = 100.0*(strtchFactor-1.0);
	std::cout << "X expansion pct. = " << strtchPct << std::endl;
	box->setSize(syssize);
	box->stretchX(strtchFactor);
	guessExpX(gel,stretchStep*nsx);
	gel->compute(true,true,false);
	output.printEnergies(gel,expXAffineFileName,strtchFactor-1.0);
	solver->solve(&model);
	gel->compute(true,true,false);
	output.printEnergies(gel,expXFileName,strtchFactor-1.0);
	//output.printSolverData(gel,solver,expXSolverFN);
	if(writeStates) {
	  fname2print.clear();
	  char expXname[128];
	  sprintf(expXname,"Run%d/expX=%f",runNum,strtchFactor);
	  fname2print.assign(expXname);
	  //gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}

	char affExpXname[128];
	sprintf(affExpXname,"Run%d/afftest-expX=%f.dat",runNum,strtchFactor);
        std::string affExpXFN(affExpXname);
        gel->storeGel(affExpXFN);
      }
    }

    if(expYtest) {
      std::cout << "Stretching gel in Y direction." << std::endl;
      for(InitPositionMap::iterator ipit=ipmap.begin(); ipit!=ipmap.end(); ipit++) {
	(ipit->first)->setPoint(ipit->second);
      }
      
      box->setShear(0.0);
      box->setSize(syssize);
      
      std::string expYFileName = nematicFileName + "-ExpY.dat";
      std::string expYAffineFileName = nematicFileName + "-affine-ExpY.dat";
      output.printParamHeader(expYFileName,paramstring.str());
      output.printParamHeader(expYAffineFileName,paramstring.str());
      //std::string fieldString = "StrY\tE_tot\tE_cl\tE_bend\tE_strc";
      std::string fieldString = "StrY\tE_tot\tE_gel\tE_bend\tE_strc\tE_cl\tE_mot\tE_pnch";
      output.printFieldLabels(expYFileName,fieldString);
      output.printFieldLabels(expYAffineFileName,fieldString);
      
      stretchStep = (stretchEnd-1.0)/nStretchSteps;
      for(int nsy=0; nsy<=nStretchSteps; nsy++) {
	double strtchFactor = 1.0+(stretchStep*nsy);
	double strtchPct = 100.0*(strtchFactor-1.0);
	std::cout << "Y expansion pct. = " << strtchPct << std::endl;
	box->setSize(syssize);
	box->stretchY(strtchFactor);
	guessExpY(gel,stretchStep*nsy);
	gel->compute(true,true,false);
	output.printEnergies(gel,expYAffineFileName,strtchFactor-1.0);
	solver->solve(&model);
	gel->compute(true,true,false);
	output.printEnergies(gel,expYFileName,strtchFactor-1.0);
	//output.printSolverData(gel,solver,expYSolverFN);
	if(writeStates) {
	  fname2print.clear();
	  char expYname[128];
	  sprintf(expYname,"Run%d/expY=%f",runNum,strtchFactor);
	  fname2print.assign(expYname);
	  //gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}

	char affExpYname[128];
	sprintf(affExpYname,"Run%d/afftest-expY=%f.dat",runNum,strtchFactor);
        std::string affExpYFN(affExpYname);
        gel->storeGel(affExpYFN);
      }
    }
    
    std::cout << "Finished testing gel." << std::endl;
    
  }
}
