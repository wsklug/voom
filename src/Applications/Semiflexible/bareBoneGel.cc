#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <getopt.h>
#include <ctime>
#include "Node.h"
#include "Gel.h"
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
#include "GelOutput_2.h"
#include "SemiflexibleInput.h"
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
  
  SemiflexibleInput inp(argv[1]);
 
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
  
  
  double actuall_c = -1.0; // actual value of l_c //
  double actualnemOP = -1.0;
  
  std::string storageFileName = inp.getString("storage file name");
  std::string storageDirectoryName = inp.getString("storage directory name");
  
  std::string fName = storageDirectoryName + storageFileName;
  
  int curGelNum = inp.getInt("gel current number");
  
  SemiflexibleGel<2> * gel;
  
  SemiflexibleGel<2>::DefNodeContainer nodes;
  PeriodicBox * box;
  GelOutput<2> output;

  // make periodic box //
  box = new LeesEdwards(inp.getReal("Wx"),inp.getReal("Wy"),0.0);

  gel = new SemiflexibleGel<2>(nodes,box,&inp);

  
  // write gel data to file //
  if(!inp.getBool("adaptiveMeshing")) {
    gel->storeGel(fName);
  }
  fName.clear();

  if(inp.getBool("relaxPrestress")) {
    gel->removePrestress();
    std::cout << "Reset spring rest lengths/stiffnesses to relax prestress." << std::endl;
  }

  
  char clinkdistfile[128];
  sprintf(clinkdistfile,"Run%d/CrosslinkSepDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,inp.getReal("L"),inp.getReal("L/l_c"),inp.getReal("dL"),inp.getReal("nematic order parameter"),curGelNum);
  fName.assign(clinkdistfile);
  output.printCrosslinkData(gel,fName);
  
  
  fName.clear();
  char fillendistfile[128];
  sprintf(fillendistfile,"Run%d/FilLengthDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,inp.getReal("L"),inp.getReal("L/l_c"),inp.getReal("dL"),inp.getReal("nematic order parameter"),curGelNum);
  fName.assign(fillendistfile);
  output.printFilLengthData(gel,fName);
  
  fName.clear();
  char nematicdistfile[128];
  sprintf(nematicdistfile,"Run%d/NematicDistro-L=%f-l_C=%f-dL=%f-S=%f-gelnum=%d.dat",runNum,inp.getReal("L"),inp.getReal("L/l_c"),inp.getReal("dL"),inp.getReal("nematic order parameter"),curGelNum);
  fName.assign(nematicdistfile);
  output.printNematicData(gel,fName);
  
  fName.clear();
  
  if(inp.getBool("ShearXTest")) {
  
    actuall_c = gel->getMeanCLsep();
    //actuall_c2 = gel2->getMeanCLsep();
    actualnemOP = gel->getNematicOP();
    double actuallamb = actuall_c*pow(actuall_c/inp.getReal("l_B"),1.0/3.0);

    char nemFileName[256];
    
    sprintf(nemFileName,"Run%d/nematicgel-L=%f-l_C=%f-dL=%f-S=%f-Wx=%f-Wy=%f-gelnum=%d",runNum,inp.getReal("L"),inp.getReal("L/l_c"),inp.getReal("dL"),inp.getReal("nematic order parameter"),curGelNum);
    
    
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
    if(inp.getString("solverType").find("LBFGSB")!=string::npos) {
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


    if(inp.getString("solverType").find("LBFGSB")!=string::npos) {
      solver = new Lbfgsb(model.dof(),m,factr,pgtol,iprint,maxiter);
    }

    else {
      solver = new Lbfgs(m,gtol,iprint,maxiter);
    }

    InitPositionMap ipmap;
    for(SemiflexibleGel<2>::DefNodeIterator n=nodes.begin(); n!=nodes.end(); n++) {
      ipmap[*n] = (*n)->point();
    }



    if(inp.getBool("zeroReturn")) {
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
    std::ostringstream paramstring;
    paramstring << "Parameters: "
		<< "L/l_C = " << inp.getReal("L")/actuall_c
		<< ", L = " << inp.getReal("L")
		<< ", S = " << actualnemOP
		<< ", l_B = " << inp.getReal("l_B")
		<< ", Wx = " << inp.getReal("Wx")
		<< ", Wy = " << inp.getReal("Wy");
    
    if(inp.getBool("checkConsist")) {
      double perturbScale = inp.getReal("minLength")/100.0;
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

    if(inp.getBool("ShearXTest")) {

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
      int nShrSteps = (int)(abs(inp.getReal("shrEnd")-inp.getReal("shrStart"))/inp.getReal("shrStep"));
      for(int sx=0; sx<=inp.getReal("shrStep"); sx++) {
	//gel->turnOffPinches();
	std::cout << "Shear X = " << inp.getReal("shrStart") + inp.getReal("shrStep")*sx << std::endl;
	box->setShearX(inp.getReal("shrStart") + inp.getReal("shrStep")*sx);
	guessAffineShearX(gel,inp.getReal("shrStart") + inp.getReal("shrStep")*sx);
	gel->compute(true,true,false);
	output.printEnergies(gel,shearXAffineFileName,inp.getReal("shrStep")*sx);
	
	solver->solve(&model);
      
	gel->compute(true,true,false);
	output.printEnergies(gel,shearXFileName,inp.getReal("shrStart")+(inp.getReal("shrStep")*sx));
	char affXname[128];
	sprintf(affXname,"Run%d/afftest-shearX=%f.dat",runNum,inp.getReal("shrStart")+(inp.getReal("shrStep")*sx));
	std::string affXFN(affXname);
	gel->storeGel(affXFN);
	//output.printSolverData(gel,solver,shearXSolverFN);
	if(inp.getBool("writeStates")) {
	  fname2print.clear();
	  char shearXname[128];
	  sprintf(shearXname,"Run%d/shearX=%f",runNum,inp.getReal("shrStart")+(inp.getReal("shrStep")*sx));
	  fname2print.assign(shearXname);
	  //gel->print(fname2print);
	  output.printGel(gel,fname2print);
	}
        
      }
    }
    std::cout << "Finished testing gel." << std::endl;
    
  }
}
