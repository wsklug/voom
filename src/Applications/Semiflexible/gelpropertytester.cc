#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <dirent.h>
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
//#include <process.h>

using namespace tvmet;
using namespace std;
using namespace voom;

typedef SemiflexibleGel<2> Gel;
typedef std::map< std::string, std::string > ParamMap;
typedef ParamMap::iterator PMIter;

typedef std::pair<double,double> doublePair;
typedef std::pair<double,doublePair> doublePairWErrors;
typedef std::map<doublePair,doublePair> dataSet;
typedef std::multimap<double,double> dataSetNoErrors;

typedef std::vector< doublePair > doublePairContainer;
typedef std::vector< doublePairWErrors > doublePairWErrorsContainer;

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

int main(int argc, char* argv[]) {

  int verbose = 1;

  char* paramFileName = argv[1];

  std::ifstream inFile(paramFileName);
  if(!(inFile.good())) {
    std::cerr << "Error: input file name does not exist!" << std::endl;
    exit(1);
  }

  char* gelFile = argv[2];
  char* unshearedgelFile;
  if(argc==5) {
    unshearedgelFile = argv[4];
  }

  double shear = atof(argv[3]);

  std::string pfn(paramFileName);
  int runnumStart = pfn.find_first_of("0123456789");
  int runnumEnd = pfn.find_first_not_of("0123456789",runnumStart);
  int runNum = atoi(pfn.substr(runnumStart,runnumEnd-runnumStart).data());

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
  double F_max = -1.0; // maximum force allowed in nonlinear rods (entropic elasticity) //
  int fitOrder = -1; // order of fitting function to use for entropic springs //
  double filDens = -1.0; // filament density in microns^-2 //
  int nNodesPerFilament = -1; // # of nodes per filament //
  double nodesPerCL = 2.5; // avg # of nodes between crosslinks //
  double nematicOP = 0.0;
  double actualnemOP = -1.0;
  std::string nemPDFparam;
  bool prestress = false;
  std::string nemPDFtype = "Gaussian";

  double k_max = -1.0; // maximum stiffness allowed in nonlinear rods (entropic elasticity) //
  double entropic_k = -1.0;
  double entropic_Lp = -1.0;

  int nGels2avg = 1;
  int curGelNum = -1;

  char gelD[16];
  sprintf(gelD,"Run%d/",runNum);

  std::string gelDirectory;
  gelDirectory.assign(gelD);
  std::string gelFileName;
  gelFileName.assign(gelFile);
  std::string bondType = "Spring";

  std::string polydisp = "none";
  double actualL;
  double minLength = 0.0;

  double motorForce = 0.0;

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
      if(parName.find("k_max")!=string::npos) {
	k_max = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("k_max",parValStr));
	bondType = "EntropicSpring";
      }
      else if(parName.find("Entropic_k")!=string::npos) {
	entropic_k = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("Entropic_k",parValStr));
      }
      else if(parName.find("Entropic_Lp")!=string::npos) {
	entropic_Lp = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("Entropic_Lp",parValStr));
      }
      else if(parName.find("kT")!=string::npos) {
	kT = atof(parValStr.data());
	pm.insert(pair< std::string, std::string>("kT",parValStr));
      }
      else if(parName.find("dt")!=string::npos) {
	dt = atof(parValStr.data());
	pm.insert(pair< std::string, std::string>("time step",parValStr));
      }
      else if(parName.find("l_B")!=string::npos) {
	l_B = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("l_B",parValStr));
      }
      else if(parName.find("kcl")!=string::npos) {
	kcl = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("crosslink stiffness",parValStr));
      }
      else if(parName.find("Max_Motor_force")!=string::npos) {
	motorForce = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("motor force",parValStr));
      }
      else if(parName.find("L_p")!=string::npos) L_p = atof(parValStr.data());
      else if(parName.find("kC")!=string::npos) kC = atof(parValStr.data());
      else if(parName.find("visc")!=string::npos) {
	visc = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("viscosity",parValStr));
      }
      else if(parName.find("MinSegLength")!=string::npos) {
        minLength = atof(parValStr.data());
      }
      else if(parName.find("L")!=string::npos) {
	if(L<0.0) L = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("L",parValStr));
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

  if(kC > 0.0 && l_B > 0.0) {
    std::ostringstream tmpStr;
    tmpStr << setprecision(16) << kC;
    pm.insert(pair< std::string, std::string>("bending modulus",tmpStr.str()));
    //gelFileName.insert(0,gelDirectory);

    gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,true,pm,minLength);
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
    
   
  }
  else {
    std::cerr << "Error: input file must have either persistence length L_p or bending modulus kC." << std::endl;
    exit(1);
  }
  
  // now prompt the user to perform various tests //

  std::cout << "Do you want to compute energy as a function of filament angle?\n(1) yes\n(2) no\n\n: ";
  int doAngleEnergy;
  std::cin >> doAngleEnergy;
  if(doAngleEnergy==1) {
    dataSet angEnerData;
    box->setShear(shear);
    angEnerData = gel->getAngularEnergyDistro();
    box->setShear(0.0);
    std::cout << "Please enter the name of a file in which to put the angular energy distribution data\n: ";
    char angFileName[256];
    std::cin >> angFileName;
    std::ofstream angFile(angFileName);
    angFile << "#tht\tE\tE_bnd\tE_strch\n";
    for(dataSet::iterator dsi=angEnerData.begin(); dsi!=angEnerData.end(); dsi++) {
      angFile << dsi->first.first << "\t" << dsi->first.second << "\t" << dsi->second.first << "\t" << dsi->second.second << "\n";
    }
    angFile.close();
  }

  std::cout << "Do you want to compute theoretical energy versus shear for independent buckling?\n(1) yes\n(2) no\n\n: ";
  int doBuckEnergy;
  std::cin >> doBuckEnergy;
  if(doBuckEnergy==1) {
    double startShear = .001;
    for(int ish=1;ish<=5;ish++) {
      std::cout << ish*startShear << "\t" << gel->computeBucklingEnergy(ish*startShear,kC,kC/sqr(l_B)) << std::endl;
    }
  }

  // COMMENTED OUT DUE TO VTK LIBRARY ISSUES ON MAC OS 10.6
  
//   std::cout << "Do you want to compute (theoretical) energy versus shear for independent buckling?\n(1) yes\n(2) no\n\n: ";
//   int doAffinityTest;
//   std::cin >> doAffinityTest;
//   if(doAffinityTest==1) {
//     std::cout << "Please input the mean crosslink distance for this gel: ";
//     double cldist;
//     std::cin >> cldist;
//     std::cout << "Please input a step size to use for averaging node displacements: ";
//     double stepSize;
//     std::cin >> stepSize;
//     std::cout << "Please input a maximum box size to use for averaging node displacements: ";
//     double maxLength;
//     std::cin >> maxLength;
//     if(cldist < minLength) cldist = minLength;
    
//     doublePairContainer affData;
//     doublePairContainer affData2;
//     box->setShear(0.0);
//     affData = gel->affineMeasurement(cldist,stepSize,maxLength,shear,"strain");
//     affData2 = gel->affineMeasurement(cldist,stepSize,maxLength,shear,"rotation");
//     std::cout << "Please enter the name of a file in which to put the affine measurement data: ";
//     char affFileName[256];
//     std::cin >> affFileName;
//     std::ofstream affFile(affFileName);
//     affFile << "#scale\tstrn\trot\n";
//     doublePairContainer::iterator dpi=affData.begin();
//     doublePairContainer::iterator dpi2=affData2.begin();
//     while(dpi!=affData.end() && dpi2!= affData2.end()) {
//       affFile << dpi->first << "\t" << dpi->second << "\t" << dpi2->second << "\n";
//       dpi++;
//       dpi2++;
//     }
//     affFile.close();
//   }


  // COMMENTED OUT DUE TO VTK LIBRARY ISSUES ON MAC OS 10.6

//   std::cout << "Do you want to compute correlation functions using the triangulated values of the rotation?\n(1) yes\n(2) no\n\n: ";
//   int doAffinityCorrTest;
//   std::cin >> doAffinityCorrTest;
//   if(doAffinityCorrTest==1) {
//     std::cout << "Please input a length scale on which to average: ";
//     double boxsize;
//     std::cin >> boxsize;
//     std::cout << "Please input a maximum length on which to compute correlations: ";
//     double maxLength;
//     std::cin >> maxLength;
    
//     std::cout << "Please input a the maximum length at which filaments are considered short: ";
//     double maxFL;
//     std::cin >> maxFL;
// //     std::cout << "Please input a file name in which to store the correlation function: ";
// //     char affCorrFN[256];
// //     std::cin >> affCorrFN;
// //     std::string acfn(affCorrFN);
    
//     box->setShear(shear);
//     gel->compute(true,true,false);
//     box->setShear(0.0);

//     gel->computeNonaffinityLengthDensityCorrelation(boxsize,maxLength,maxFL,shear);
//   }

  std::cout << "Do you want to compute affinity as a function of length scale (Head/Levine measure)?\n(1) yes\n(2) no\n\n: ";
  int doAffinityTest2;
  std::cin >> doAffinityTest2;
  if(doAffinityTest2==1) {
    std::cout << "Please input the mean crosslink distance for this gel: ";
    double cldist;
    std::cin >> cldist;
    std::cout << "Please input a step size to use for averaging node displacements: ";
    double stepSize;
    std::cin >> stepSize;
    std::cout << "Please input a maximum box size to use for averaging node displacements: ";
    double maxLength;
    std::cin >> maxLength;
    if(cldist < minLength) cldist = minLength;
    std::cout << "Please input a maximum size of filaments to consider: ";
    double maxFilSize;
    std::cin >> maxFilSize;
    std::cout << "Do you want to use pairs of points or points and interpolated points?\n(1) point pairs\n(2) points/interpolated points\n\n: ";
    int interp;
    std::cin >> interp;
    bool doAngleTest = false;
    std::cout << "Do you want to keep track of the nonaffinity as a function of angle?\n(1) yes\n(2) no\n\n: ";
    int doa;
    std::cin >> doa;
    if(doa==1) doAngleTest = true;

    doublePairWErrorsContainer affData3;
    box->setShear(0.0);
    
    if(interp==1) {
      affData3 = gel->affineMeasurementHeadLevine(cldist,stepSize,maxLength,shear,maxFilSize,doAngleTest);
    }
    else {
      affData3 = gel->affineMeasurementHeadLevineInterpolated(cldist,stepSize,maxLength,shear,maxFilSize);      
    }
    std::cout << "Please enter the name of a file in which to put the Head/Levine affine measurement data: ";
    char affFileName2[256];
    std::cin >> affFileName2;
    std::ofstream affFile2(affFileName2);
    affFile2 << "#r/L\tdth^2\terr\n";
    doublePairWErrorsContainer::iterator dpi3=affData3.begin();
    while(dpi3!=affData3.end()) {
      affFile2 << (dpi3->first)/L << "\t" << dpi3->second.first << "\t" << dpi3->second.second << "\n";
      dpi3++;
    }
    affFile2.close();
  }

  std::cout << "Do you want to compute affinity as a function of position (Head/Levine measure) for a fixed length scale?\n(1) yes\n(2) no\n\n: ";
  int doAffinityTest3;
  std::cin >> doAffinityTest3;
  if(doAffinityTest3==1) {
    std::cout << "Please input a pair separation distance: ";
    double pairDist;
    std::cin >> pairDist;
    std::cout << "Please input a box size: ";
    double boxSize;
    std::cin >> boxSize;
    if(boxSize <= pairDist) {
      std::cout << "Sorry, box size must be equal to or greater than pair separation distance.  Resetting box size to " << 2.0*pairDist << "." << std::endl;
      boxSize = 2.0*pairDist;
    }

    tvmet::Vector<double,2> boxDims;
    boxDims[0] = boxSize;
    boxDims[1] = boxSize;

    std::cout << "Please input a maximum size of filaments to consider: ";
    double maxFilSize;
    std::cin >> maxFilSize;

    box->setShear(0.0);
   
    std::cout << "Please enter a file name in which to put the nonaffinity data (leave out the vtk extension):";
    char affFileName3[256];
    std::cin >> affFileName3;
    std::string aFN3(affFileName3);
    
    std::cout << "Would you also like to compute the nonaffinity correlation function? ";
    int nacf;
    std::cin >> nacf;
    bool doNACorr = false;
    if(nacf==1) doNACorr = true;
    
    gel->affineBoxesMeasurement(boxDims,pairDist,shear,maxFilSize,aFN3,doNACorr);
  }


  std::cout << "Do you want to compute nonaffinity/energy/length density cross-correlations?\n(1) yes\n(2) no\n\n: ";
  int doCrossCorr;
  std::cin >> doCrossCorr;
  if(doCrossCorr==1) {
    std::cout << "Please enter a length scale on which to mesh the system: ";
    double len;
    std::cin >> len;
    std::cout << "Please enter a file name in which to put the cross-correlation data: ";
    char crossCorrFN[256];
    std::cin >> crossCorrFN;
    std::string crossCorrFile(crossCorrFN);
    box->setShear(shear);
    gel->compute(true,true,false);
    box->setShear(0.0);
    gel->computeCrossCorrelations(len,shear,crossCorrFile);
  }

  std::cout << "Do you want to make a buckling order parameter (MacKintosh/Conti measure) map?\n(1) yes\n(2) no\n\n: ";
  int doBuckOPMap;
  std::cin >> doBuckOPMap;
  if(doBuckOPMap==1) {
    std::cout << "Please enter a length scale on which to mesh the system: ";
    double len;
    std::cin >> len;
    bool goodVal = true;
    std::map<double,std::string> shearGelMap;
    shearGelMap.insert(pair<double,std::string>(0.0,gelFileName));
    while(goodVal) {
      std::cout << "Please enter the shear exactly as it appears in the next sheared gel file's name, or simply type DONE to indicate all shear files are exhausted:\n";
      char shearedFN[256];
      std::cin >> shearedFN;
      std::string shearFNandShear(shearedFN);
      if(shearFNandShear.find("DONE")!=string::npos) goodVal = false;
      else {
	double newShear = atof(shearedFN);
	int spaceEq = gelFileName.find_last_of("=");
	int spaceDot = gelFileName.find_last_of(".");
	std::string shearFN(gelFileName);
	shearFN.replace(spaceEq+1,spaceDot-spaceEq-1,shearFNandShear);
	shearGelMap.insert(pair<double,std::string>(newShear,shearFN));
      }
    }

    tvmet::Vector<double,2> gridSize;
    gridSize[0] = len;
    gridSize[1] = len;

    std::cout << "Do you want to compute correlations (this will take more time)?\n(1) yes\n(2) no\n\n: ";
    int doCorr;
    std::cin >> doCorr;
    bool doCorrs = false;
    if(doCorr==1) doCorrs = true;

    gel->buckleOPCalc(len,shearGelMap,doCorrs);

    box->setShear(0.0);
  }

  std::cout << "Do you want to make a buckling map?\n(1) yes\n(2) no\n\n: ";
  int doBuckMap;
  std::cin >> doBuckMap;
  if(doBuckMap==1) {
    std::cout << "Please enter a length scale on which to mesh the system: ";
    double len;
    std::cin >> len;
    std::cout << "At what fraction of bending energy should compressed segments be considered buckled? ";
    double bfrac;
    std::cin >> bfrac;
    std::cout << "Do you want to compute the buckling correlation function, too? ";
    int dobc;
    std::cin >> dobc;
    bool compCF = true;
    if(dobc!=1) compCF = false;
    std::cout << "Please enter a file name in which to store the buckling map (leave out the vtk extension): ";
    char crossCorrFN[256];
    std::cin >> crossCorrFN;
    std::string crossCorrFile(crossCorrFN);
    crossCorrFile += ".vtk";
    box->setShear(shear);
    gel->compute(true,true,false);
    box->setShear(0.0);

    tvmet::Vector<double,2> gridSize;
    gridSize[0] = len;
    gridSize[1] = len;

    gel->computeBucklingMap(gridSize,shear,bfrac,crossCorrFile,compCF);

    box->setShear(0.0);
  }

  std::cout << "Do you want to compute the buckling strain map and correlation?\n(1) yes\n(2) no\n\n: ";
  int doBuckStrainMap;
  std::cin >> doBuckStrainMap;
  if(doBuckStrainMap==1) {
    std::cout << "Please enter a length scale on which to mesh the system: ";
    double len;
    std::cin >> len;
    std::cout << "At what fraction of bending energy should compressed segments be considered buckled? ";
    double bfrac;
    std::cin >> bfrac;
    bool goodVal = true;
    std::map<double,std::string> shearGelMap;
    //shearGelMap.insert(pair<double,std::string>(0.0,gelFileName));
    while(goodVal) {
      std::cout << "Please enter the shear exactly as it appears in the next sheared gel file's name, or simply type DONE to indicate all shear files are exhausted:\n";
      char shearedFN[256];
      std::cin >> shearedFN;
      std::string shearFNandShear(shearedFN);
      if(shearFNandShear.find("DONE")!=string::npos) goodVal = false;
      else {
	double newShear = atof(shearedFN);
	int spaceEq = gelFileName.find_last_of("=");
	int spaceDot = gelFileName.find_last_of(".");
	std::string shearFN(gelFileName);
	shearFN.replace(spaceEq+1,spaceDot-spaceEq-1,shearFNandShear);
	shearGelMap.insert(pair<double,std::string>(newShear,shearFN));
      }
    }

    tvmet::Vector<double,2> gridSize;
    gridSize[0] = len;
    gridSize[1] = len;

    std::cout << "Do you want to compute correlations (this will take more time)?\n(1) yes\n(2) no\n\n: ";
    int doCorr;
    std::cin >> doCorr;
    bool doCorrs = false;
    if(doCorr==1) doCorrs = true;

    gel->cooperativeBuckleMeasure(len,bfrac,l_B,shearGelMap,doCorrs);

    box->setShear(0.0);
  }

  std::cout << "Do you want to compute the fraction of energy in filaments parallel and perpendicular to the shear?\n(1) yes\n(2) no\n\n: ";
  int doParPerpTest;
  std::cin >> doParPerpTest;
  if(doParPerpTest==1) {
    box->setShear(shear);
    gel->compute(true,false,false);
    double fpar = gel->parallelenergy()/(gel->energy());
    double fperp = gel->perpenergy()/(gel->energy());
    std::cout << "Gel energy = " << gel->energy() << ", f_parallel = " << fpar << ", f_perp = " << fperp << "." << std::endl;
    box->setShear(0.0);
  }

  std::cout << "Do you want to print out the angles of all filaments?\n(1) yes\n(2) no\n\n: ";
  int doAngTest;
  std::cin >> doAngTest;
  if(doAngTest==1) {
    std::cout << "Please enter the name of a file in which to store the angles: ";
    char angF[256];
    std::cin >> angF;
    std::string aF(angF);
    gel->printAngles(aF);
  }

  std::cout << "Do you want to check for large bends?\n(1) yes\n(2) no\n\n: ";
  int doBendTest;
  std::cin >> doBendTest;
  if(doBendTest==1) {
    box->setShear(shear);
    gel->compute(true,false,false);
    gel->printBigBends();
    box->setShear(0.0);
  }

 std::cout << "Do you want to compute the density/energy correlation function?\n(1) yes\n(2) no\n\n: ";
  int doDensTest;
  std::cin >> doDensTest;
  if(doDensTest==1) {
    box->setShear(shear);
    gel->compute(true,true,false);
    box->setShear(0.0);
    std::cout << "Please enter a length scale to mesh the system: ";
    double gridSize;
    std::cin >> gridSize;
    std::multimap< double, vector<double> > densData;
    densData = gel->getDensityEnergyDistro(gridSize);
    std::cout << "Please enter the name of a file in which to put the density/energy distribution data\n: ";
    char densFileName[256];
    std::cin >> densFileName;
    std::ofstream densFile(densFileName);
    densFile << "#dens\tEtot\tEbend\tEstrch\n";
    for(multimap< double, vector<double> >::iterator mi=densData.begin(); mi!=densData.end(); mi++) {
      densFile << mi->first << "\t" << mi->second[0] << "\t" << mi->second[1] << "\t" << mi->second[2] << "\n";
    }
    densFile.close();
  }

 std::cout << "Do you want to compute the energy correlation function?\n(1) yes\n(2) no\n\n: ";
  int doEnCorrTest;
  std::cin >> doEnCorrTest;
  if(doEnCorrTest==1) {
    box->setShear(shear);
    gel->compute(true,true,false);
    box->setShear(0.0);
    std::cout << "Please enter a length scale to mesh the system: ";
    double gridSize;
    std::cin >> gridSize;
    tvmet::Vector<double,2> gridDims;
    gridDims[0] = gridSize;
    gridDims[1] = gridSize;
    std::cout << "Please enter a maximum separation to consider: ";
    double maxl;
    std::cin >> maxl;
    std::vector< pair<double, double> > enCorrDat;
    enCorrDat = gel->energyCorrelationFunction(gridSize,maxl);
    std::cout << "Please enter the name of a file in which to put the energy correlation data\n: ";
    char enCorrFileName[256];
    std::cin >> enCorrFileName;
    std::ofstream enCorrFile(enCorrFileName);
    enCorrFile << "#r\tcorr\n";
    for(vector< pair<double,double> >::iterator corri=enCorrDat.begin(); corri!=enCorrDat.end(); corri++) {
      enCorrFile << corri->first << "\t" << corri->second << "\n";
    }
    enCorrFile.close();
  }

  std::cout << "Do you want to compute the energies of affine deformations?\n(1) yes\n(2) no\n\n: ";
  int doAffDefTest;
  std::cin >> doAffDefTest;
  if(doAffDefTest==1) {
    std::vector< pair<double,double> > shearXAffData;
    box->setShear(0.0);
    gel->compute(true,true,false);
    std::pair<double,double> datPair = pair<double,double>(0.0,gel->energy());
    shearXAffData.push_back(datPair);
    std::cout << "Please enter a shear step size: ";
    double shearStepSize;
    std::cin >> shearStepSize;
    std::cout << "Please enter the number of steps to take (excluding shear=0.0): ";
    int nShearSteps;
    std::cin >> nShearSteps;
    for(int nsx=1; nsx<=nShearSteps; nsx++) {
      guessAffineShearX(gel,nsx*shearStepSize);
      box->setShearX(nsx*shearStepSize);
      gel->compute(true,true,false);
      datPair = pair<double,double>(nsx*shearStepSize,gel->energy());
      double bende = gel->bendingenergy();
      if(abs(bende) > 1.0e-8) std::cout << "Sanity check: bending energy = " << bende << " at shear = " << nsx*shearStepSize << "." << std::endl;
      shearXAffData.push_back(datPair);
      box->setShearX(0.0);
    }
    std::cout << "Please enter the name of a file in which to put the affine energy data\n: ";
    char affEnFileName[256];
    std::cin >> affEnFileName;
    std::ofstream affEnFile(affEnFileName);
    affEnFile << "#shrX\ten\n";
    for(vector< pair<double,double> >::iterator pi=shearXAffData.begin(); pi!=shearXAffData.end(); pi++) {
      affEnFile << pi->first << "\t" << pi->second << std::endl;
    }
    affEnFile.close();
  }

  std::cout << "Do you want to test the nonlinear springs?\n(1) yes\n(2) no\n\n: ";
  int doNLST;
  std::cin >> doNLST;
  if(doNLST == 1) {
    std::vector<int> dofim(2);
    dofim[0] = 0;
    dofim[1] = 1;
    int id = 0;
    tvmet::Vector<double,2> p1,p2;
    p1[0] = 0.0;
    p1[1] = 0.0;
    p2[0] = 0.1;
    p2[1] = 0.0;
    double L0 = norm2(p2-p1);
    id++;
    dofim[0] += 2;
    dofim[1] += 2;
    BrownianNode<2> * node1 = new BrownianNode<2>(id,dofim,p1,p1);
    BrownianNode<2> * node2 = new BrownianNode<2>(id,dofim,p2,p2);

    std::cout << "L0 = " << L0 << std::endl
	      << "Lp = ";
    
    double Lp;
    std::cin >> Lp;
    
    std::cout << "entrop_k = ";

    double entropk;
    std::cin >> entropk;

    std::cout << "mu = ";
    double mu;
    std::cin >> mu;

    EntropicSpring<2> * newESpring = new EntropicSpring<2>(node1,node2,entropk,Lp,L0,mu);

    std::cout << "File in which to store spring data: ";

    char fn[256];
    std::cin >> fn;

    std::string fname(fn);

    newESpring->printOutVals(fname);

    newESpring->checkConsistency();

    delete newESpring;
  }

  std::cout << "Do you want to compute the energy of the gel?\n(1) yes\n(2) no\n\n: ";
  int doEnergyCalc;
  std::cin >> doEnergyCalc;
  if(doEnergyCalc==1) {
    box->setShear(shear);
    gel->compute(true,true,false);
    std::cout << shear << "\t" << gel->energy() << "\t" << gel->filenergy() << "\t" << gel->bendingenergy() << "\t" << gel->stretchingenergy()<< "\t" << gel->crosslinkenergy() << "\t" << gel->motorenergy() << "\t" << gel->pinchenergy() << std::endl;
    box->setShear(0.0);
  }

  std::cout << "Do you want to minimize the energy in an old gel?\n(1) yes\n(2) no\n\n: ";
  int doEnergyMinCalc;
  std::cin >> doEnergyMinCalc;
  if(doEnergyMinCalc==1) {
    box->setShear(0.0);
    gel->compute(true,true,false);
    std::cout << "Sanity check: gel energy at zero shear = " << gel->energy() << std::endl;
    Model::NodeContainer modelnodes;
    modelnodes.reserve( nodes.size() );
    for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      modelnodes.push_back( *n );
    }
    Model model(modelnodes);
    model.pushBackBody( gel );
    
    int m=7;
    double factr=1.0e7;
    double pgtol=1.0e-5;
    double gtol = 1.0e-5;
    int iprint = 100;
    int maxiter = 1500000;
    Solver * solver;

//     ifstream lbfgsbinp("lbfgsb.inp");
//     lbfgsbinp >> iprint >> factr >> pgtol >> m ;
//     if(verbose) {
//       std::cout << "Input iprint: " << iprint << std::endl
// 		<< "Input factr: " << factr << std::endl
// 		<< "Input pgtol: " << pgtol << std::endl
// 		<< "Input m: " << m << std::endl;
//     }
    solver = new Lbfgsb(model.dof(),m,factr,pgtol,iprint,maxiter);

    std::cout << "Enter a shear value: ";
    double shearin;
    std::cin >> shearin;

    box->setShear(shearin);
    guessAffineShearX(gel,shearin);
    solver->solve(&model);
    
    gel->compute(true,true,false);
    std::cout << shearin << "\t" << gel->energy() << "\t" << gel->filenergy() << "\t" << gel->bendingenergy() << "\t" << gel->stretchingenergy()<< "\t" << gel->crosslinkenergy() << "\t" << gel->motorenergy() << "\t" << gel->pinchenergy() << std::endl;
    box->setShear(0.0);
  }

//   std::cout << "Do you want to run a consistency check on entropic springs?\n(1) yes\n(2) no\n: ";
//   int doConsCheck;
//   std::cin >> doConsCheck;
//   if(doConsCheck && bondType == "EntropicSpring") {
//     SemiflexibleGel<2>::Filament * Fil;
//     int nFils = gel->filaments().size();
//     for(int numfil=0; numfil<nFils; numfil++) {
//       Fil = gel->filament(numfil);
//       int nESprings = Fil->bonds.size();
//       for(int numsp=0; numsp<nESprings; numsp++) {
// 	bool consist = Fil->bonds[numsp]->checkConsistency();
// 	if(!consist) std::cout << "Filament " << numfil << ", bond " << numsp << " is inconsistent!" << std::endl;
//       }
//     }
//   }

  std::cout << "Do you want to make a solution movie/see when nodes converged?\n(1) yes\n(2) no\n: ";
  int doSolMovie;
  std::cin >> doSolMovie;
  if(doSolMovie==1) {
    char runDirName[256];
    sprintf(runDirName,"Run%d",runNum);
    DIR* runDir = opendir(runDirName);
    std::map<int,std::string> stateFiles;
    if(runDir!=0) {
      struct dirent* newFile = readdir(runDir);
      while(newFile!=0) {
	char* stateFileName = newFile->d_name;
	std::string stateFileNameStr(stateFileName);
	if(stateFileNameStr.find("lbfgsbiter")!=string::npos) {
	  int iterStart = stateFileNameStr.find_first_of("0123456789");
	  int iterEnd = stateFileNameStr.find_first_not_of("0123456789",iterStart);
	  int iterNumber = atoi(stateFileNameStr.substr(iterStart,iterEnd-iterStart).data());
	  stateFiles.insert(pair<int,std::string>(iterNumber,stateFileNameStr));
	}
	newFile = readdir(runDir);
      }

      std::cout << "Read in data on " << stateFiles.size() << " solution snapshots." << std::endl;

      std::vector< tvmet::Vector<double,2> > finalPositions;
      int nGelFils = gel->filaments().size();
      for(int ngf=0; ngf<nGelFils; ngf++) {
	int nGelNodes = gel->filament(ngf)->nodes.size();
	for(int ngn=0; ngn<nGelNodes; ngn++) {
	  finalPositions.push_back(gel->filament(ngf)->nodes[ngn]->point());
	}
      }
      
      std::map<int,std::string>::iterator sfiter = stateFiles.begin();
      for(sfiter; sfiter!=stateFiles.end(); sfiter++) { 
	gelFileName = sfiter->second;
	//gelFileName.insert(0,gelDirectory);
	SemiflexibleGel<2>::DefNodeContainer newNodes;
	SemiflexibleGel<2>* newGel = new SemiflexibleGel<2>(gelFileName,newNodes,bondType,true,pm,minLength);

	newGel->setBox(box);
	newGel->makeFinalPosMap(finalPositions);

	char gelStoreName[256];
	sprintf(gelStoreName,"lbfgsb-iter%d",sfiter->first);
	std::string gelSN(gelStoreName);
	box->setShear(shear);
	newGel->compute(true,true,false);
	std::cout << "Gel energy = " << newGel->energy() << std::endl;
	//gel->print(gelStoreName);
	output.printGel(newGel,gelStoreName);
	box->setShear(0.0);
	delete newGel;	
      }
      
    }

    else {
      std::cout << "Error: could not open directory " << runDirName << " with solution snapshots." << std::endl;
    }

  }

  std::cout << "Do you want to write out segment connectivity data?\n(1) yes\n(2) no\n: ";
  int doSegConn;
  std::cin >> doSegConn;
  if(doSegConn==1) {
    std::cout << "Please enter a file name in which to store the segment connectivity data: ";
    char gelConnName[256];
    std::cin >> gelConnName;
    std::string gelCN(gelConnName);

    gel->printConnectionData(gelCN);
  }

  std::cout << "Do you want to write out a list of segments that have stiffened by at least some amount?\n(1) yes\n(2) no\n: ";
  int doSegWrite;
  std::cin >> doSegWrite;
  if(doSegWrite==1) {
    double minStiff;
    std::cout << "By what minimum factor should filaments have stiffened to be considered for inclusion? ";
    std::cin >> minStiff;

    char stiffSegFile[256];
    std::cout << "Enter the name of a file in which to store the stiff segments: ";
    std::cin >> stiffSegFile;
    std::string stiffSF(stiffSegFile);
    
    box->setShear(shear);
    gel->compute(true,true,false);

    gel->printStiffenedSegments(stiffSF,minStiff);
  }


  std::cout << "Do you want to store the state of the gel to a vtk file for later Paraview visualization?\n(1) yes\n(2) no\n: ";
  int doStore;
  std::cin >> doStore;
  if(doStore==1) {
    std::cout << "Please enter a file name in which to store the gel's state (do not provide an extension, please): ";
    char gelStoreName[256];
    std::cin >> gelStoreName;
    std::string gelSN(gelStoreName);
    box->setShear(shear);
    gel->compute(true,true,false);
    std::cout << "Gel energy = " << gel->energy() << std::endl;
    //gel->print(gelStoreName);
    output.printGel(gel,gelStoreName);
    box->setShear(0.0);
  }

  std::cout << "Do you want to create an abaqus input file for this gel?\n(1) yes\n(2) no\n: ";
  int doAba;
  std::cin >> doAba;
  if(doAba==1) {
    std::cout << "Please enter a file name in which to store the gel's state for abaqus (do not provide an extension, please): ";
    char gelStoreName[256];
    std::cin >> gelStoreName;
    std::string gelSN(gelStoreName);
    gelSN += ".inp";
    box->setShear(0.0);
    std::cout << "Please enter the value of mu to use: ";
    char muchar[256];
    std::cin >> muchar;
    std::cout << "Please enter the value of l_B to use: ";
    char lBchar[256];
    std::cin >> lBchar;
    double mu4ab = atof(muchar);
    double lB4ab = atof(lBchar);
    std::cout << "Do you want to create nodes at the edges or use a 'ghost' node to enforce proper displacements?\n(1) edge pin\n(2) ghost node\n:";
    int method;
    std::cin >> method;
    double minL;
    if(method==1) {
      std::cout << "Please enter a minimum segment length: ";
      std::cin >> minL;
    }
    else minL = 0.0;
    gel->printforAbaqus(gelSN,mu4ab,lB4ab,minL,method);
  }
  
}
