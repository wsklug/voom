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
//#include "Lbfgs.h"
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
      if(parName.find("kT")!=string::npos) {
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
    gelFileName.insert(0,gelDirectory);

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
  
  std::cout << "Do you want to compute (theoretial) energy versus shear for independent buckling?\n(1) yes\n(2) no\n\n: ";
  int doAffinityTest;
  std::cin >> doAffinityTest;
  if(doAffinityTest==1) {
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
    
    doublePairContainer affData;
    doublePairContainer affData2;
    box->setShear(0.0);
    affData = gel->affineMeasurement(cldist,stepSize,maxLength,shear,"strain");
    affData2 = gel->affineMeasurement(cldist,stepSize,maxLength,shear,"rotation");
    std::cout << "Please enter the name of a file in which to put the affine measurement data: ";
    char affFileName[256];
    std::cin >> affFileName;
    std::ofstream affFile(affFileName);
    affFile << "#scale\tstrn\trot\n";
    doublePairContainer::iterator dpi=affData.begin();
    doublePairContainer::iterator dpi2=affData2.begin();
    while(dpi!=affData.end() && dpi2!= affData2.end()) {
      affFile << dpi->first << "\t" << dpi->second << "\t" << dpi2->second << "\n";
      dpi++;
      dpi2++;
    }
    affFile.close();
  }

  std::cout << "Do you want to compute correlation functions using the triangulated values of the rotation?\n(1) yes\n(2) no\n\n: ";
  int doAffinityCorrTest;
  std::cin >> doAffinityCorrTest;
  if(doAffinityCorrTest==1) {
    std::cout << "Please input a length scale on which to average: ";
    double boxsize;
    std::cin >> boxsize;
    std::cout << "Please input a maximum length on which to compute correlations: ";
    double maxLength;
    std::cin >> maxLength;
    
    std::cout << "Please input a the maximum length at which filaments are considered short: ";
    double maxFL;
    std::cin >> maxFL;
//     std::cout << "Please input a file name in which to store the correlation function: ";
//     char affCorrFN[256];
//     std::cin >> affCorrFN;
//     std::string acfn(affCorrFN);
    
    box->setShear(shear);
    gel->compute(true,true,false);
    box->setShear(0.0);

    gel->computeNonaffinityLengthDensityCorrelation(boxsize,maxLength,maxFL,shear);
  }

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


  std::cout << "Do you want to compute cross-correlations?\n(1) yes\n(2) no\n\n: ";
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
    gel->print(gelStoreName);
    box->setShear(0.0);
  }
  
}
