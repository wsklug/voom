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
#include "BrownianRod.h"
#include "BrownianDynamics.h"
#include <random/uniform.h>
#include <random/normal.h>
//#include "Dirichlet.h"
#include "PeriodicBox.h"
#include "LeesEdwards.h"
#include "GelOutput.h"
//#include "Grid.h"
//#include "TwoBodyPotential.h"
//#include <process.h>

#define USEOPENMP 0
#if defined(_OPENMP)
#include <omp.h>
#undef USEOPENMP
#define USEOPENMP 1
#endif

using namespace tvmet;
using namespace std;
using namespace voom;

typedef SemiflexibleGel<2> Gel;
typedef SemiflexibleGel<2>::Filament Fil;
typedef SemiflexibleGel<2>::DefNode BNode;
typedef SemiflexibleGel<2>::DefNodeContainer BNodeContainer;
typedef std::map<std::string,std::string> ParamMap;
typedef ParamMap::iterator ParamIterator;
typedef std::map< Fil*,std::set<Fil*> > PairContainer;
typedef PairContainer::iterator PairIterator;

int main(int argc, char* argv[]) {

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
  sprintf(newDirName,"mkdir -p RodLiquidRun%d",runNum);
  system(newDirName);

  double kT = -1.0; // temperature in pN-microns //
  double visc = -1.0; // viscosity //
  tvmet::Vector<double,2> syssize; // system size in microns //
  double kBond = -1.0; // spring constant for rods in pN/micron //
  double dt = -1.0; // time step in seconds //
  double L = -1.0; // length of filaments in microns //
  double filDens = -1.0; // filament density in microns^-2 //
  
  enum TwoBodyPotential::FunctionType ft;
  ft = TwoBodyPotential::gaussian;
  double maxPPsep = -1.0;
  double chargeSep = -1.0;
  double energyScale = -1.0;
  double lengthScale = -1.0;

  double maxTime = -1.0;
  double printTime = -1.0;

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
      else if(parName.find("visc")!=string::npos) {
	visc = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("viscosity",parValStr));
      }
      else if(parName.find("dt")!=string::npos) {
	dt = atof(parValStr.data()); 
	pm.insert(pair< std::string, std::string >("time step",parValStr));
      }
      else if(parName.find("Wx")!=string::npos) syssize[0] = atof(parValStr.data());
      else if(parName.find("Wy")!=string::npos) syssize[1] = atof(parValStr.data());

      else if(parName.find("FilDens")!=string::npos) filDens = atof(parValStr.data());

      else if(parName.find("kBond")!=string::npos) {
	kBond = atof(parValStr.data());
	pm.insert(pair<std::string, std::string>("filament spring constant",parValStr));
      }

      else if(parName.find("Charge Sep.")!=string::npos) {
	chargeSep = atof(parValStr.data());
	pm.insert(pair<std::string, std::string>("separation of \"charges\"",parValStr));
      }
      else if(parName.find("Max. PP Sep.")!=string::npos) {
	maxPPsep = atof(parValStr.data());
      }
      else if(parName.find("EnergyScale")!=string::npos) {
	energyScale = atof(parValStr.data());
	pm.insert(pair<std::string,std::string>("energy scale for two body potential",parValStr));
      }
      else if(parName.find("LengthScale")!=string::npos) {
	lengthScale = atof(parValStr.data());
	pm.insert(pair<std::string,std::string>("length scale for two body potential",parValStr));
      }
      else if(parName.find("MaxTime")!=string::npos) {
	maxTime = atof(parValStr.data());
      }
      else if(parName.find("PrintTime")!=string::npos) {
	printTime = atof(parValStr.data());
      }
      else if(parName.find("PotentialType")!=string::npos) {
	if(parValStr.find("Exp")!=string::npos) ft = TwoBodyPotential::expon;
      }
      else if(parName.find("L")!=string::npos) {
	L = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("filament length",parValStr));
      }      
    }
    curString.clear();
    parName.clear();
    parValStr.clear();
  }
  inFile.close();

  if(maxPPsep < 0.0) maxPPsep = 2.0*L;

  if(chargeSep < 0.0) chargeSep = .005;

  int nChrgs = (int)(L/chargeSep);
  chargeSep = L/nChrgs;
  nChrgs++;

  double vol = 1.0;
  for(int i=0;i<2;i++) {
    syssize[i] *= L;
    vol *= syssize[i];
  }
  
  std::ostringstream ssstream;
  ssstream << "(" << syssize[0] << ", " << syssize[1] << ")";
  pm.insert(pair<std::string,std::string>("system size",ssstream.str()));

  int nFils = (int)(filDens*vol);
  filDens = nFils/vol;
  ssstream.str("");
  ssstream.flush();
  ssstream << filDens;
  pm.insert(pair<std::string, std::string>("filament density",ssstream.str()));
  ssstream.str("");
  ssstream.flush();
  ssstream << nFils;
  pm.insert(pair<std::string, std::string>("number of filaments",ssstream.str()));

  Gel * gel = new SemiflexibleGel<2>();
  PeriodicBox * box = new PeriodicBox(syssize[0],syssize[1]);
  gel->setBox(box);
  
  Vector2D gridSpaces;
  gridSpaces[0] = maxPPsep/(2.0*sqrt(2.0));
  gridSpaces[1] = maxPPsep/(2.0*sqrt(2.0));
  Grid<Fil,Fil,2> * grid = new Grid<Fil,Fil,2>(box,gridSpaces,&Fil::point,true);
  gridSpaces = grid->gridSpace();
  maxPPsep = gridSpaces[0]*4.0*sqrt(2.0);
  ssstream.str("");
  ssstream.flush();
  ssstream << maxPPsep/2.0;
  pm.insert(pair<std::string,std::string>("maximum sep. for two body force calculation",ssstream.str()));
  ssstream.str("");
  ssstream.flush();
  ssstream << gridSpaces[0];
  pm.insert(pair<std::string, std::string>("grid spacing",ssstream.str()));
  
  TwoBodyPotential * tbp = new TwoBodyPotential(box,ft,maxPPsep,energyScale/sqr(nChrgs),L,lengthScale,nChrgs);
//   Vector2D R;
//   double thet;
//   int it=0;
//   while(it<10) {
//     std::cout << "Testing two body potential.\nEnter Rx: ";
//     std::cin >> R[0];
//     std::cout << "Enter Ry:";
//     std::cin >> R[1];
//     std::cout << "Enter theta:";
//     std::cin >> thet;
//     std::cout << "Energy = " << tbp->integrateEnergy(R,thet) << std::endl;
//     std::cout << "Energy = " << tbp->integrateForce(R,thet) << std::endl;
//     std::cout << "Energy = " << tbp->integrateTorque(R,thet) << std::endl;
//     it++;
//   }
 
  ranlib::Uniform<double> rnguni;
  rnguni.seed((unsigned int)time(0));
  int id = 0;
  NodeBase::DofIndexMap idx(2);
  BNodeContainer nodes;
  for(int i=0; i<nFils; i++) {
    BNodeContainer bns;
    for(int j=0; j<2; j++) {
      idx[j] = 2*id + j;
    }
    rnguni.seed((unsigned int)time(0)+i);
    Vector2D pos;
    pos[0] = rnguni.random()*syssize[0]*L;
    pos[1] = rnguni.random()*syssize[1]*L;
    BNode * newNode = new BrownianNode<2>(id,idx,pos,pos);
    newNode->setId(id);
    id++;
    bns.push_back(newNode);
    nodes.push_back(newNode);
    double newAng = 2.0*M_PI*rnguni.random();
    Vector2D endpos;
    endpos[0] = pos[0] + cos(newAng)*L;
    endpos[1] = pos[1] + sin(newAng)*L;
    for(int j=0; j<2; j++) {
      idx[j] = 2*id + j;
    }
    newNode = new BrownianNode<2>(id,idx,endpos,endpos);
    newNode->setId(id);
    id++;
    bns.push_back(newNode);
    nodes.push_back(newNode);
    gel->addFilament(bns,kBond,1.0,visc,kT,dt);
  }

  std::cout << "Adding filaments to grid..." << std::endl;
  std::vector<Fil*> fils = gel->filaments();
  grid->addElems(fils);
  
//   for(int testcount=0; testcount<20; testcount++) {
//     int ri = (int)(rnguni.random()*nFils);
//     Fil* f = gel->filament(ri);
//     PairContainer & pairs = grid->getNeighbors();
//     int nn = pairs[f].size();
//     int rn = (int)(rnguni.random()*nn);
//     set<Fil*>::iterator fi = pairs[f].begin();
//     for(int j=0;j<rn;j++) fi++;
//     Fil* f2 = *fi;
//     tbp->checkTable(f->nodes,f2->nodes);
//   }

  int maxIters = (int)(maxTime/dt);
  int printStep = (int)(printTime/dt);

  BrownianDynamics* bd = new BrownianDynamics(nodes,-1);
  gel->setGrid(grid);
  gel->addTwoBodyPotential(tbp);
  bd->pushBackBody(gel);

  GelOutput<2> outp;

  char infofl[256];
  sprintf(infofl,"RodLiquidRun%d/runinfo.dat",runNum);
  std::ofstream infofile(infofl);
  for(ParamIterator pi=pm.begin(); pi!=pm.end(); pi++) {
    infofile << pi->first << "\t" << pi->second << std::endl;
  }
  infofile.close();
  
  char datfl[256];
  sprintf(datfl,"RodLiquidRun%d/nemdata.dat",runNum);
  std::ofstream dfile(datfl);
  dfile << "#t\tS_tot\tang" << std::endl;
  dfile.close();

  for(int ct=0; ct<maxIters/printStep; ct++) {
    bd->computeAndAssemble(true,false,false);
    std::cout << "Brownian Dynamics: " << ct*printStep << " time steps taken; time = " << ct*printStep*dt << ", energy = " << gel->energy() << std::endl;
    char fname[256];
    sprintf(fname,"RodLiquidRun%d/gelliq-%d",runNum,ct);
    std::string fnm(fname);
    outp(gel,fnm);
    gel->computeNematicOP();
    Vector2D & director = gel->getNemDirector();
    double nemang = atan2(director[1],director[0]);
    dfile.open(datfl,ios_base::app);
    dfile << ct*printStep*dt << "\t" << gel->getNematicOP() << "\t" << nemang << std::endl;
    dfile.close();
    char nemfl[256];
    sprintf(nemfl,"RodLiquidRun%d/nemcorrdata%d.dat",runNum,ct);
    std::vector< std::pair<double,double> > nemcorrdata = gel->computeNemCorrelations(.5*L,.5*L,.1*L);
    std::ofstream nemf(nemfl);
    nemf << "#r\tcorr" << std::endl;
    for(std::vector< std::pair<double,double> >::iterator di=nemcorrdata.begin(); di!=nemcorrdata.end(); di++) {
      nemf << di->first << "\t" << di->second << std::endl;
    }
    nemf.close();
    
    bd->run(printStep,dt);
  }
  
  char messg[256];
  sprintf(messg,"Finished with rod liquid simulation; output found in directory RodLiquidRun%d\n",runNum);
  std::string endmessg(messg);
  std::cout << endmessg << std::endl;
}
