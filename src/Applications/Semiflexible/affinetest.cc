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
  char * gelStoreFN = new char[128];
  double L = atof(pmap["L"].data());
  double l_c = L/(atof(pmap["L/l_c"].data()));
  double kcl = atof(pmap["crosslink stiffness"].data());
  double Wx = atof(pmap["Wx"].data());
  double Wy = atof(pmap["Wy"].data());
  double dL = atof(pmap["dL"].data());
  std::string gelStoreName;
  
  if(gelNum == -1){
    sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-gelnum=1.gelsave",L,l_c,dL,Wx,Wy,kcl);
    gelStoreName.assign(gelStoreFN);
    gelStoreName.insert(0,gelLibDir);
    gelNum = getCurGelNum(gelStoreName);
    gelStoreName.clear();
    delete gelStoreFN;
    gelStoreFN = new char[128];
  }
  
  sprintf(gelStoreFN,"gel-L=%f-l_C=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-gelnum=%d.gelsave",L,l_c,dL,Wx,Wy,kcl,gelNum);
  gelStoreName.assign(gelStoreFN);
  gelStoreName.insert(0,gelLibDir);
  return gelStoreName;
}

// void setParamsfromFileName(std::string fileName, ParamMap & pmap) {
  
// }

double getNematicParam(double nemOP) {
  // takes nematic order parameter and returns proper parameters for prob. dist. //
  double a0 = .0613966;
  double a1 = 3.49418;
  double a2 = -2.63328;
  double num = a0 + a1*nemOP + a2*sqr(nemOP);
  return num/(1.0-nemOP);
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
  double lambda = -1.0; // l_c*(l_c/l_B)^z //
  double F_max = -1.0; // maximum force allowed in nonlinear rods (entropic elasticity) //
  int fitOrder = -1; // order of fitting function to use for entropic springs //
  double filDens = -1.0; // filament density in microns^-2 //
  int nNodesPerFilament = -1; // # of nodes per filament //
  double nodesPerCL = 2.5; // avg # of nodes between crosslinks //
  double nematic = 0.0;
  bool prestress = false;
  bool cutoffends = false;
  
  bool affineTest = false;
  std::string affmethod;
  double affshear = -1.0;
  double affbinspace = 0.2;
  double affbintol = affbinspace*(1.0e-3);
  double affmin = affbinspace;
  double affmax = -1.0;
  int nGels2avg = 1;
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
      else if(parName.find("Affine_test") != string::npos) {
	if(atof(parValStr.data()) >= 0.5) affineTest = true;
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
      else if(parName.find("l_B")!=string::npos) {
	l_B = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("l_B",parValStr));
      }
      else if(parName.find("Lambda")!=string::npos) lambda = atof(parValStr.data());
      else if(parName.find("fit")!=string::npos) {
	fitOrder = atoi(parValStr.data());
	pm.insert(pair< std::string, std::string >("fit order",parValStr));	
      }
      else if(parName.find("Affine_method")!=string::npos) affmethod.assign(parValStr);
      else if(parName.find("Affine_shear")!=string::npos) affshear = atof(parValStr.data());
      else if(parName.find("Affine_gels")!=string::npos) nGels2avg = atoi(parValStr.data());
      else if(parName.find("Affine_space")!=string::npos) affbinspace = atof(parValStr.data());
      else if(parName.find("Affine_tol")!=string::npos) affbintol = atof(parValStr.data());
      else if(parName.find("Affine_min")!=string::npos) affmin = atof(parValStr.data());
      else if(parName.find("Affine_max")!=string::npos) affmax = atof(parValStr.data());
      else if(parName.find("Nematic")!=string::npos) {
	nematic = atof(parValStr.data());
	if(nematic > 1.0e-6) {
	  std::ostringstream nem;
	  nem << getNematicParam(nematic);
	  pm.insert(pair< std::string, std::string >("nematic order parameter",nem.str()));
	}
	else {
	  pm.insert(pair< std::string, std::string >("nematic order parameter",parValStr));
	}
      }
      else if(parName.find("Prestress")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) pm.insert(pair< std::string, std::string >("prestress","false"));
	else pm.insert(pair< std::string, std::string >("prestress","true"));
      }
      else if(parName.find("Retrieve")!=string::npos) {
	if(atof(parValStr.data()) < 0.5) retrieveGel = false;
	else retrieveGel = true;
      }
      else if(parName.find("L")!=string::npos) {
	L = atof(parValStr.data());
	pm.insert(pair< std::string, std::string >("L",parValStr));
      }
      else if(parName.find("CutOffEnds")!=string::npos) {
	if(atof(parValStr.data()) >= 0.5) cutoffends = true; 
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
      pm["L"] = Lstr;
      L = atof(Lstr.data());
      std::string l_Cstr;
      l_Cstr.assign(gelFileName,l_Cpos+4,dLpos-5-l_Cpos);
      tmpratio = L/atof(l_Cstr.data());
      std::ostringstream tm;
      tm << tmpratio;
      pm["L/l_c"] = tm.str();
      std::string dLstr;
      dLstr.assign(gelFileName,dLpos+3,Wxpos-4-dLpos);
      pm["dL"] = dLstr;
      dL = atof(dLstr.data());
      std::string kclStr;
      kclStr.assign(gelFileName,kclpos+4,gelNumpos-5-kclpos);
      kcl = atof(kclStr.data());
      pm["crosslink stiffness"] = kclStr;
      std::string GNstr;
      GNstr.assign(gelFileName,gelNumpos+7,endPos-gelNumpos-7);
      curGelNum = atoi(GNstr.data());
      gelFileName.insert(0,gelDirectory);
      if(lambda  > 0.0) { // if lambda is specified, override l_B //
	l_B = pow(L/tmpratio,4)/pow(lambda,3);
	std::ostringstream tmplB;
	tmplB << l_B;
	pm["l_B"] = tmplB.str();
      }
      else {
	lambda = (tmpratio/L)*pow(tmpratio/(L*l_B),1/3);
      }
      std::ostringstream lambstr;
      lambstr << lambda;
      pm["lambda"] = lambstr.str();
      gel = new SemiflexibleGel<2>(gelFileName,nodes,bondType,cutoffends,pm);
      box = gel->box();
      syssize = box->size();
      std::cout << "Retrieved and set up a gel with the following properties:" << std::endl;
      std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
      std::cout << "Bond type: " << bondType << std::endl;
      for(PMIter pmi=pm.begin(); pmi!=pm.end(); pmi++) {
	std::cout << pmi->first << ": " << pmi->second << std::endl;
      }
      std::ostringstream ss0;
      ss0 << syssize[0];
      std::ostringstream ss1;
      ss1 << syssize[1];
      pm["Wx"] = ss0.str();
      pm["Wy"] = ss1.str();
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

    if(lambda  > 0.0) { // if lambda is specified, override l_B //
      l_B = pow(L/tmpratio,4)/pow(lambda,3);
      std::ostringstream tmplB;
      tmplB << l_B;
      pm["l_B"] = tmplB.str();
      kBond = kAngle/sqr(l_B);
      std::ostringstream tmpStr;
      tmpStr << kBond;
      pm["bond stiffness"] = tmpStr.str();
    }
    else {
      lambda = (tmpratio/L)*pow(tmpratio/(L*l_B),1.0/3.0);
    }
    std::ostringstream lambstr;
    lambstr << lambda;
    pm["lambda"] = lambstr.str();
    
    std::ostringstream tmpSx;
    tmpSx << syssize[0];
    pm.insert(pair< std::string, std::string >("Wx",tmpSx.str()));
    std::ostringstream tmpSy;
    tmpSy << syssize[1];
    pm.insert(pair< std::string, std::string >("Wy",tmpSy.str()));
    std::ostringstream tmpLl_c;
    tmpLl_c << tmpratio;
    pm.insert(pair< std::string, std::string >("L/l_c",tmpLl_c.str()));
    std::ostringstream tmpdL;
    tmpdL << dL;
    pm.insert(pair< std::string, std::string >("dL",tmpdL.str()));

    std::cout << "Constructing a gel with the following properties:" << std::endl;
    std::cout << "System size: " << syssize[0] << ", " << syssize[1] << std::endl;
    std::cout << "# nodes/filament: " << nNodesPerFilament << std::endl;
    std::cout << "Bond type: " << bondType << std::endl;
    for(PMIter pmi=pm.begin(); pmi!=pm.end(); pmi++) {
      std::cout << pmi->first << ": " << pmi->second << std::endl;
    }

    // make periodic box //
    box = new LeesEdwards(syssize[0],syssize[1],0.0);
    
    // create body //
    gel = new SemiflexibleGel<2>(nodes,box,filDens,nNodesPerFilament,dL,bondType,cutoffends,pm);
    //gel->addPinch(6.0,false,nodes,kBond,kAngle,visc,kT,dt,kcl);
    
    // write gel data to file //
    std::string fName = getGelFileName(gelDirectory,pm);
    curGelNum = getCurGelNum(fName);
    gel->storeGel(fName);

    fName.clear();
    char clinkdistfile[128];
    sprintf(clinkdistfile,"CrosslinkSepDistro-L=%f-l_C=%f-dL=%f-gelnum=%d.dat",L,L/tmpratio,dL,curGelNum);
    fName.assign(clinkdistfile);
    fName.insert(0,gelDirectory);
    std::cout << "Printing out crosslink distribution." << std::endl;
    output.printCrosslinkData(gel,fName);
  }

  if(affineTest) {
    actuall_c = gel->getMeanCLsep();
    double actuallamb = actuall_c*pow(actuall_c/l_B,1.0/3.0);
    double shear = 0.0;
    if(affmethod.find("Langer")!=string::npos) {

    }
    else {
      // first, create arrays to hold data //
      if(affmax<0.0) {
	affmax = (syssize[0]+syssize[1])/(4.0*L);
      }
      int nPts = (int)((affmax-affmin)/affbinspace);
      DataSet data(nPts);
      for(int pt=0; pt<nPts; pt++) {
	data[pt].rL = affmin + affbinspace*pt;
	data[pt].nDatPts = 0;
	data[pt].mean = 0.0;
	data[pt].stddev = 0.0;
      }
      double rMax = data[nPts-1].rL + .5*affbinspace;
      // now set up model and solve with no shear //
      Model::NodeContainer modelnodes;
      modelnodes.reserve( nodes.size() );
      std::cout << "# of nodes in model = " << nodes.size() << std::endl;
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
    
      output( gel, "relaxed-0" );
      
      solver.solve( &model );  

      // now collect data on node pair angles //
      NodePairList nodePairs[nPts];
//       for(int j=0; j<nPts; j++) {
// 	NodePairList tmpNPL;
// 	nodePairs.push_back(tmpNPL);
//       }
      // iterate through all filament pairs and add relevant node pair data to data array //
      int nFils = gel->filaments().size();
      int f1,f2;
      for(f1=0; f1<nFils; f1++) {
	Gel::Filament * fil1 = gel->filament(f1);
	for(Gel::ConstDefNodeIterator ni1=(fil1->nodes).begin(); ni1!= (fil1->nodes).end(); ni1++) {
	  tvmet::Vector<double,2> n1vec;
	  n1vec = (*ni1)->point();
	  if(box->inside(n1vec)) {
	    for(Gel::ConstDefNodeIterator ni1p=ni1+1; ni1p!=(fil1->nodes).end(); ni1p++) {
	      //if(!(gel->isSlave(*ni1p))) {
	      tvmet::Vector<double,2> n1pvec;
	      n1pvec = (*ni1p)->point();
	      if(box->inside(n1pvec)) {
		tvmet::Vector<double,2> initDiff;
		initDiff = n1vec - n1pvec;
		double r = norm2(initDiff)/L;
		if(r<rMax && r>=affmin-affbintol && fmod(r-affmin+affbintol,affbinspace) <= 2.0*affbintol) {
		  // determine proper bin for data point //
		  double tmpInd = (r-affmin+1.5*affbintol)/affbinspace;
		  int indx = (int)(tmpInd);
		  
		  // compute affine value of Deltatheta for pair of points //
		  double thetaInit = atan2(initDiff[1],initDiff[0]);
		  NodePair newNP;
		  newNP.node1 = *ni1;
		  newNP.node2 = *ni1p;
		  newNP.initAngle = thetaInit;
		  nodePairs[indx].push_back(newNP);
		}
	      }
	    }
	  
	    for(f2=f1+1; f2<nFils; f2++) {
	      Gel::Filament * fil2 = gel->filament(f2);
	      for(Gel::ConstDefNodeIterator ni2=(fil2->nodes).begin(); ni2!=(fil2->nodes).end(); ni2++) {
		//if(!(gel->isSlave(*ni2))) {
		if(!(gel->areLinked(*ni1,*ni2))) {
		  tvmet::Vector<double,2> n2vec;
		  n2vec = (*ni2)->point();
		  if(box->inside(n2vec)) {
		    tvmet::Vector<double,2> initDiff;
		    initDiff = n1vec - n2vec;
		    double r = norm2(initDiff)/L;
		    if(r<rMax && r>=affmin-affbintol && fmod(r-affmin+affbintol,affbinspace) <= 2.0*affbintol) {
		      // determine proper bin for data point //
		      double tmpInd = (r-affmin+1.5*affbintol)/affbinspace;
		      int indx = (int)(tmpInd);
		      
		      // compute affine value of Deltatheta for pair of points //
		      double thetaInit = atan2(initDiff[1],initDiff[0]);
		      NodePair newNP;
		      newNP.node1 = *ni1;
		      newNP.node2 = *ni2;
		      newNP.initAngle = thetaInit;
		      nodePairs[indx].push_back(newNP);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      // now shear box and solve model again //
      box->setShear(affshear);
      shear = affshear;
      solver.solve(&model);

      // collect data on final angles and add differences to data list, then clear out angle data and node pair lists //
      for(int k=0; k<nPts; k++) {
	NodePairList::iterator npl;
	for(npl=nodePairs[k].begin(); npl!=nodePairs[k].end(); npl++) {
	  Gel::DefNode* n1 = npl->node1;
	  Gel::DefNode* n2 = npl->node2;
	  tvmet::Vector<double,2> npt1;
	  tvmet::Vector<double,2> npt2;
	  npt1 = n1->point();
	  npt2 = n2->point();
	  if(box->inside(npt1) && box->inside(npt2)) {
	    double thetaInit = npl->initAngle;
	    tvmet::Vector<double,2> initDiff;
	    initDiff[0] = cos(thetaInit);
	    initDiff[1] = sin(thetaInit);
	    double thetaAff = atan2(initDiff[1],initDiff[0]+affshear*initDiff[1]);
	    tvmet::Vector<double,2> finalDiff;
	    finalDiff = npt1 - npt2;
	    double thetaFinal = atan2(finalDiff[1],finalDiff[0]);
	    // 	  if(thetaFinal-thetaInit < 0.0) dThetaDiff = thetaFinal - atan2(initDiff[1],initDiff[0]+affshear*initDiff[1]);
	    // 	  else dThetaDiff = -2.0*M_PI + thetaFinal - atan2(initDiff[1],initDiff[0]+affshear*initDiff[1]);
	    
	    double dThetaDiff = (fabs((thetaFinal - thetaAff)) < 2*M_PI - fabs((thetaFinal - thetaAff))) ? thetaFinal - thetaAff : 2*M_PI - fabs(thetaFinal - thetaAff); 
	    double newDatPt = sqr(dThetaDiff);
	    // 	  if(newDatPt > 1.0) {
	    // 	    std::cerr << "Large angle deviation: r/L = " << data[k].rL << ", DThetaAffine = " << atan2(initDiff[1],initDiff[0]+affshear*initDiff[1])-thetaInit << ", DTheta = " << thetaFinal-thetaInit << std::endl;
	    // 	    std::cerr << "Node 1 located at (" << npt1[0] << "," << npt1[1] << "); Node 2 located at (" << npt2[0] << "," << npt2[1] << ")." << std::endl << std::endl;
	    // 	  }
	    data[k].mean += newDatPt; // update running total //
	    data[k].stddev += sqr(newDatPt); // update running total of squared value //
	    data[k].nDatPts++; // update total number of data points for this distance //
	  } 
	}
	nodePairs[k].clear();
      }

      for(int nGel=1; nGel<nGels2avg; nGel++) {
	if(retrieveGel) {
	  // if first gel was retrieved, check for another gel; if it is there, read it in, and if not, create a new one //
	}
	else {
	  // if first gel was created anew, create another gel //
	}
	// repeat steps from above on nth gel //
      }

      // calculate mean, etc. from data //
      for(DataSet::iterator di=data.begin(); di!=data.end(); di++) {
	int nDatPts = di->nDatPts;
	double mean = (di->mean)/nDatPts;
	mean /= sqr(affshear);
	double stddev = (di->stddev/(nDatPts*sqr(sqr(affshear))))-sqr(mean);
	double fct = ((double)(nDatPts))/((double)(nDatPts-1));
	stddev *= fct;
	stddev = sqrt(stddev);
	di->mean = mean;
	di->stddev = stddev/sqrt(nDatPts);
      }

      // output data to file //
      char affFileName[128];
      sprintf(affFileName,"affinemeasure-L=%f-l_C=%f-l_B=%f-lambda=%f-dL=%f-Wx=%f-Wy=%f-kcl=%f-gelnum=%d.dat",L,L/tmpratio,l_B,lambda,dL,syssize[0],syssize[1],kcl,curGelNum);
      std::ostringstream paramstring;
      paramstring << "#Parameters: "
		  << "Affine measure = " << "Head/Levine"
		  << ", shear = " << affshear
		  << ", actual L/l_c = " << L/actuall_c
		  << ", actual lambda = " << actuallamb
		  << ", r/L tol = " << affbintol << std::endl;
      
      std::string fieldHeaders = "#r/L\t<Dth^2>\tstd. dev.\tN_pts";
      ofstream affFile(affFileName);
      affFile << paramstring.str() << fieldHeaders << std::endl;
      for(DataSet::iterator di=data.begin(); di!=data.end(); di++) {
	affFile << di->rL << "\t" << di->mean << "\t" << di->stddev << "\t" << di->nDatPts << std::endl;
      }
      affFile.close();
    } 
  }
}
