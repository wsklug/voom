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
#include "BrownianRod.h"
#include <random/uniform.h>
#include <random/normal.h>
//#include "Dirichlet.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "GelOutput.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

typedef std::map< std::string, std::string > ParamMap;
typedef ParamMap::iterator PMIter;

double getFun(double nu) {
  double num = nu - 1.0 + exp(-nu);
  double denom = 1.0 + exp(-nu) - (2.0/nu)*(1.0 - exp(-nu));
  return num/denom;
}

double getFilDens(double l_C, double L) {
  double f = L/l_C;
  double nulow = 0.0;
  double nuhigh = f;
  double curf = getFun((nulow+nuhigh)/2.0);
  while(abs(curf-f) > 1.0e-6) {
    if(curf >= f) nuhigh = (nuhigh+nulow)/2.0;
    else nulow = (nuhigh+nulow)/2.0;
    curf = getFun((nulow+nuhigh)/2.0);
  }
  return M_PI*((nulow+nuhigh)/4.0)/sqr(L);
}

int getFilSegs(double avgCL, double nodesPerCL) {
  return (int)((avgCL+1)*nodesPerCL);
}

int main(int argc, char* argv[]) {

  int verbose = 1;

  char* paramFileName = argv[1];

  std::ifstream inFile(paramFileName);
  if(!(inFile.good())) {
    std::cerr << "Error: input file name does not exist!" << std::endl;
    exit(1);
  }

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
  double tmpratio = -1.0; // L/l_C //
  double actuall_c = -1.0; // actual value of l_c //
  double actualL = -1.0;
  double F_max = -1.0; // maximum force allowed in nonlinear rods (entropic elasticity) //
  int fitOrder = -1; // order of fitting function to use for entropic springs //
  double filDens = -1.0; // filament density in microns^-2 //
  int nNodesPerFilament = -1; // # of nodes per filament //
  double nodesPerCL = 2.5; // avg # of nodes between crosslinks //
  bool prestress = false; // whether to prestress filaments (put in wiggles) //
  bool cutOffEnds = false; // whether to cut off ends of filaments past last crosslinks //
  bool minBeforeStore = false; // whether to run a minimization before storing the gel //
  
  std::string polydisp = "None";

  bool shearTest = false;
  double shrStart = -1.0;
  double shrEnd = -1.0;
  int nShrSteps = -1;
  double shrStep = -1.0;

  int curGelNum = -1;

  std::string gelDirectory = "./";
  std::string gelFileName;
  std::string bondType = "Spring";

  bool retrieveGel = false;

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
      else if(parName.find("Shear_test") != string::npos) {
	if(atof(parValStr.data()) >= 0.5) shearTest = true;
      }
      else if(parName.find("GelDirectory") != string::npos) gelDirectory.assign(parValStr);
      else if(parName.find("kT")!=string::npos) {
	kT = atof(parValStr.data());
	pm.insert(pair< std::string, std::string>("kT",parValStr));
      }
      else if(parName.find("L_p")!=string::npos) L_p = atof(parValStr.data());
      else if(parName.find("kC")!=string::npos) kC = atof(parValStr.data());
      else if(parName.find("visc")!=string::npos) {
	visc = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("viscosity",parValStr));
      }
      // else if(parName.find("r")!=string::npos) r = parVal;
      else if(parName.find("F_max")!=string::npos) {
	F_max = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("F_max",parValStr));
	bondType = "EntropicSpring";
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
      else if(parName.find("l_B")!=string::npos) l_B = atof(parValStr.data());
      else if(parName.find("fit")!=string::npos) {
	fitOrder = atoi(parValStr.data());
	pm.insert(pair< std::string, std::string >("fit order",parValStr));	
      }
      else if(parName.find("Shear_start")!=string::npos) shrStart = atof(parValStr.data());
      else if(parName.find("Shear_final")!=string::npos) shrEnd = atof(parValStr.data());
      else if(parName.find("Shear_steps")!=string::npos) nShrSteps = atoi(parValStr.data());
      else if(parName.find("Prestress")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) pm.insert(pair< std::string, std::string >("prestress","false"));
	else pm.insert(pair< std::string, std::string >("prestress","true"));
      }
      else if(parName.find("Retrieve")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) retrieveGel = false;
	else retrieveGel = true;
      }
      else if(parName.find("CutOffEnds")!=string::npos) {
	if(atof(parValStr.data()) > .5) cutOffEnds = true;
      }
      else if(parName.find("Presave min")!=string::npos) {
	if(atof(parValStr.data()) > .5) minBeforeStore = true;
      }
      else if(parName.find("Polydispersity")!=string::npos) {
	polydisp = parValStr;
	pm.insert(pair< std::string, std::string >("polydispersity",polydisp));
      }
      else if(parName.find("Long/short")!=string::npos && polydisp.find("LongShort")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longshortratio",parValStr));
      }
      else if(parName.find("Long fraction")!=string::npos && polydisp.find("LongShort")!=string::npos) {
	pm.insert(pair< std::string, std::string >("longfraction",parValStr));
      }
      else if(parName.find("Min. length")!=string::npos && polydisp.find("Exponential")!=string::npos) {
	pm.insert(pair< std::string, std::string >("minlength",parValStr));
      }
      else if(parName.find("L_fil")!=string::npos) {
	L = atof(parValStr.data());
      } 
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
      tmpStr << kC;
      pm.insert(pair< std::string, std::string>("bending modulus",tmpStr.str()));
      std::ostringstream tmpStr2;
      tmpStr2 << l_B;
      pm.insert(pair< std::string, std::string>("l_B",tmpStr2.str()));
      // get nominal l_c and L from gel file name //
      int Lpos = gelFileName.find("L=") + 2;
      int l_Cpos = gelFileName.find("l_C");
      int dLpos = gelFileName.find("dL");
      int Wxpos = gelFileName.find("Wx");
      int kclpos = gelFileName.find("kcl");
      int gelNumpos = gelFileName.find("gelnum");
      int endPos = gelFileName.find(".gelsave");
      std::string Lstr;
      Lstr.assign(gelFileName,Lpos,l_Cpos-1-Lpos);
      L = atof(Lstr.data());
      std::string l_Cstr;
      l_Cstr.assign(gelFileName,l_Cpos+4,dLpos-5-l_Cpos);
      tmpratio = L/atof(l_Cstr.data());
      std::string dLstr;
      dLstr.assign(gelFileName,dLpos+3,Wxpos-4-dLpos);
      dL = atof(dLstr.data());
      std::string kclStr;
      kclStr.assign(gelFileName,kclpos+4,gelNumpos-5-kclpos);
      kcl = atof(kclStr.data());
      pm["crosslink stiffness"] = kclStr;
      std::string GNstr;
      GNstr.assign(gelFileName,gelNumpos+7,endPos-gelNumpos-7);
      curGelNum = atoi(GNstr.data());
      gelFileName.insert(0,gelDirectory);
      gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutOffEnds,pm);
      box = gel->box();
      syssize = box->size();
    }
    else {
      std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
      exit(1);
    }
  }

  else {
    if(tmpratio >= 0.0 && L >= 0.0) {
      filDens = getFilDens(L/tmpratio,L);
      syssize[0] *= L;
      syssize[1] *= L;
    }
    else {
      std::cerr << "Error: crosslink density information not found; please check input file " << paramFileName << " to see whether L and L/l_C are defined."  << std::endl;
      std::cerr << "L = " << L << ", L/l_C = " << tmpratio << std::endl;
      exit(1);
    }
    
    // get number of rods per filament from requirement of nodes between crosslinks, then set all other parameters //
    nNodesPerFilament = 1 + getFilSegs(tmpratio-1,nodesPerCL);
    dL = L/(nNodesPerFilament - 1);
    // Estimate kAngle = E*I/L = kT\xi_p/dL; kT=4.1pN-nm, \xi_p=10^4nm, dL~100nm
    //make kAngle 100 times larger to simulate the rotation diffusion
    
    if(kC > 0.0) {
      kAngle = kC/dL;
      std::ostringstream tmpStr;
      tmpStr << kAngle;
      pm.insert(pair< std::string, std::string>("angle stiffness",tmpStr.str()));
    }
    else {
      std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
      exit(1);
    }
    
    if(l_B > 0.0) {
      kBond = kAngle/sqr(l_B);
      std::ostringstream tmpStr;
      tmpStr << kBond;
      pm.insert(pair< std::string, std::string>("bond stiffness",tmpStr.str()));
    }
    else {
      std::cerr << "Error: input file must have a value for the bending/stretching length l_B." << std::endl;
    }
    box = new LeesEdwards(syssize[0],syssize[1],0.0);
    // create body
    std::cout << "Constructing a gel with the following properties:" << std::endl;
    std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
    std::cout << "# nodes/filament: " << nNodesPerFilament << std::endl;
    std::cout << "Bond type: " << bondType << std::endl;
    for(PMIter pmi=pm.begin(); pmi!=pm.end(); pmi++) {
      std::cout << pmi->first << ": " << pmi->second << std::endl;
    }
    
    gel = new SemiflexibleGel<2>(nodes,box,filDens,nNodesPerFilament,dL,bondType,cutOffEnds,pm);
    //gel->addPinch(6.0,false,nodes,kBond,kAngle,visc,kT,dt,kcl);

    std::ifstream testStream;
    bool fileFound = true;
    char * gelStoreFN;
    int i = 1;
    std::string fName;
    
    // find first available gel name //
    while(fileFound) {
      fName.clear();
      gelStoreFN = new char[128];
      sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-gelnum=%d.gelsave",L,L/tmpratio,dL,syssize[0],syssize[1],kcl,i);
      fName.assign(gelStoreFN);
      fName.insert(0,gelDirectory);
      testStream.open(fName.c_str());
      if(testStream.fail()) {
	fileFound = false;
	curGelNum = i;
      }
      delete gelStoreFN;
      i++;
      testStream.close();
    }

    char clinkdistfile[128];
    sprintf(clinkdistfile,"CrosslinkSepDistro-L=%f-l_C=%f-dL=%f-gelnum=%d.dat",L,L/tmpratio,dL,curGelNum);
    std::string clfName;
    clfName.assign(clinkdistfile);
    clfName.insert(0,gelDirectory);
    std::cout << "Printing out crosslink distribution." << std::endl;
    output.printCrosslinkData(gel,clfName);
    
    if(polydisp.find("Exponential")!=string::npos || polydisp.find("LongShort")!=string::npos) {
      char fillendistfile[128];
      sprintf(fillendistfile,"FilamentLengthDistro-L=%f-l_C=%f-dL=%f-gelnum=%d.dat",L,L/tmpratio,dL,curGelNum);
      std::string fName2(fillendistfile);
      fName2.insert(0,gelDirectory);
      std::cout << "Printing out filament length distribution" << std::endl;
      output.printFilamentLengthData(gel,fName2);
    }
   
    if(!minBeforeStore) {
      gel->storeGel(fName);
    }
    else {
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
      int iprint = 0;
      ifstream lbfgsbinp("lbfgsb.inp");
      lbfgsbinp >> iprint >> factr >> pgtol >> m ;
      int maxiter = 1500000;
      Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxiter );
      solver.solve(&model);
      gel->storeGel(fName);
    }
  }
  
  // check to make sure that no two filaments are multiply linked //
  // gel->checkCrosslinks();
  
  if(shearTest) {
    actuall_c = gel->getMeanCLsep();
    if(polydisp.find("Exponential")!=string::npos || polydisp.find("LongShort")!=string::npos) {
      actualL = gel->getMeanFilLen();
    }
    else actualL = L;
    
    if(shrStart >= 0.0 && shrEnd >= shrStart && nShrSteps > 0) {
      shrStep = (shrEnd-shrStart)/nShrSteps;
      std::cout << "About to begin system shear from " << shrStart << " to " << shrEnd << " in steps of size " << shrStep << std::endl;
    }
    else {
      std::cerr << "Error: input file must list starting and ending shear and # of shear steps to take." << std::endl;
    }
    
    
    //   // Estimate kBond = EA/L; EA~ 12 EI/d^2 = 12*\xi_p*kT/d^2; d~8nm;
    //   // EA ~ 2/3 * 10^4 pN
    //   double r=3.5e-3; // radius in microns
    //   double kBond = (1.0e-4)*4.0*kAngle/(r*r); // pN/micron
    
    //   // 1 Pa = 1 N/m^2 = 10^12 pN/(10^9 nm)^2 = 10^-6 pN/nm^2 = 10^-6 MPa
    //   // \eta_s = 10^-3 Pa-s = 10^-9 pN-s/nm^2
    // //   double viscosity = 1.0e-9; // MPa-s
    //   double viscosity = 0.0;//1.0; // pN-ms/micron^2
    
    //   // drag = 2\pi*\eta_s*L; L=100nm, 
    //   // drag = 6*10^-7 pN-s/nm
    
    // //   double dt = 1.0e-7; // s
    //   double dt = 2.5e-5; // ms
    
    //   double kcl = -kBond/10.0;
    
    //   double filDens = 14.0;
    
    //   double maxForce = 100.0; // pN
    
    char shrfilename[128];
    sprintf(shrfilename,"energyvsshear-L=%f-l_C=%f-l_B=%f-gelnum=%d.dat",L,L/tmpratio,l_B,curGelNum);
    
    
    //   std::cout << "Derived parameters:" << std::endl
    // 	    << "# of nodes/filament = " << nNodesPerFilament << std::endl
    // 	    << "kAngle = " << kAngle << " pN-microns" << std::endl
    // 	    << "kBond = " << kBond << " pN/micron" << std::endl
    // 	    << "Wx = " << syssize[0] << " microns" << std::endl
    // 	    << "Wy = " << syssize[1] << " microns" << std::endl;
    
    
    
    /////////////////////////////////////////////////////////////////////////////
    //
    // Create a model and a minimization solver and point them to the
    // SemiflexibleGel body.
    //
    /////////////////////////////////////////////////////////////////////////////
    
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
    int iprint = 0;
    ifstream lbfgsbinp("lbfgsb.inp");
    lbfgsbinp >> iprint >> factr >> pgtol >> m ;
    if(verbose) 
      std::cout << "Input iprint: " << iprint << std::endl
		<< "Input factr: " << factr << std::endl
		<< "Input pgtol: " << pgtol << std::endl
		<< "Input m: " << m << std::endl;
    int maxiter = 1500000;
    Lbfgsb solver(model.dof(), m, factr, pgtol, iprint, maxiter );
    
    std::string shearFileName(shrfilename);
    std::ostringstream paramstring;
    paramstring << "Parameters: "
	      << "L/l_C = " << actualL/actuall_c
	      << ", L = " << actualL
	      << ", l_B = " << l_B
	      << ", Wx = " << syssize[0]
	      << ", Wy = " << syssize[1];
    output.printParamHeader(shearFileName,paramstring.str());

    std::ostringstream energyFLstring;
    
    energyFLstring << "Shear\t" << "E_tot\t" << "E_cl\t" << "E_bend\t" << "E_strc";
    output.printFieldLabels(shearFileName,energyFLstring.str());
    
    for(int step=0; step<=nShrSteps; step++) {
      double shear=shrStart+step*shrStep;
      cout << "Step " << step << " shear = " << shear << endl;
      box->setShear(shear);
      solver.solve( &model );
      char fname[100];
      sprintf(fname,"shear-%d",step);
      output(gel, fname);
      gel->print(fname);
      output.printEnergies(gel,shearFileName,shear);      
    }
  }

  std::cout << "Completed shear test." << std::endl;
  
  return 0; 

}
