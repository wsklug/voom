// -*- C++ -*-
//----------------------------------------------------------------------
//
//                          Andrew R. Missel
//                University of California Los Angeles
//                   (C) 2010 All Rights Reserved
//
//----------------------------------------------------------------------

#if !defined(__ClusterFinder_h__)
#define __ClusterFinder_h__

#include <set>
#include <string>
#include <map>
#include <multimap.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>

class ClusterFinder {
  
 public:

  struct Site {
    bool active;
    std::vector<double> loc;
    int id;
    std::vector<int> neighbs;
  };
  

  struct Cluster {
    std::set<Site*> members;
  };

  ClusterFinder(std::string & clusterDefFile, std::vector<double> & size) { 
    _h = size;

    std::cout << "Read in data on system size: (" << size[0] << "," << size[1] << ")" << std::endl;

    std::ifstream clustFile(clusterDefFile.c_str());

    int id = 0;
    while(!clustFile.eof()) {
      bool gooddata = true;
      char curLine[256];
      clustFile.getline(curLine,256);
      std::string curString;
      curString.assign(curLine);
      int breakchar = curString.find_first_of(",");
      std::string tmpString;
      tmpString = curString.substr(0,breakchar);
      if(breakchar == std::string::npos) {
	gooddata = false;
      }
      else {
	Site* newSite = new Site();
	newSite->id = id;
	id++;
	newSite->active = true;
	
	newSite->loc.push_back(atof(tmpString.data()));
	int oldbreak = breakchar;
	breakchar = curString.find_first_of(",",oldbreak+1);
	newSite->loc.push_back(atof((curString.substr(oldbreak+1,breakchar-oldbreak-1)).data()));
	oldbreak = breakchar;
	breakchar = curString.find_first_of(",",oldbreak+1);
	while(breakchar!=std::string::npos) {
	  newSite->neighbs.push_back(atoi((curString.substr(oldbreak+1,breakchar-oldbreak-1)).data()));
	  oldbreak = breakchar;
	  breakchar = curString.find_first_of(",",oldbreak+1);
	}
	_allsites.push_back(newSite);
	//std::assert(_allsites.size() == id);
      }
      if(_allsites.size() != id) {
	std::cout << "There is something rotten in the state of Denmark..." << std::endl;
      }
    }

    clustFile.close();

    std::cout << "Read in data on " << _allsites.size() << " sites." << std::endl;
  }

  ~ClusterFinder() {

  }

  void readActivity(std::string & activityFile) {
    std::vector<Site*>::iterator si = _allsites.begin();
    for(si; si!=_allsites.end(); si++) {
      (*si)->active = false;
    }
    
    std::ifstream actFile(activityFile.c_str());

    while(!actFile.eof()) {
      char curLine[256];
      actFile.getline(curLine,256);
      std::string curString;
      curString.assign(curLine);
      if(curString.find_first_of("0123456789")!=std::string::npos) {
	_allsites[atoi(curString.data())]->active = true;
      }
    }

    actFile.close();
  }

  void mapDistance(std::vector<double> & dX) {
    for(int i=0; i<2; i++) {
      dX[i] = dX[i] - 
	( dX[i]>0.0 ? std::floor((dX[i]/_h[i])+0.5) : std::ceil((dX[i]/_h[i])-0.5))*_h[i];
    }      
    return;
  }

  void addToCluster(Site* newSite, Cluster* cluster) {
    if(newSite->active && _clustersites.find(newSite) == _clustersites.end()) {
      _clustersites.insert(newSite);
      cluster->members.insert(newSite);
      std::vector<int>::iterator ni = newSite->neighbs.begin();
      for(ni; ni!=newSite->neighbs.end(); ni++) {
	addToCluster(_allsites[*ni],cluster);
      }
    }
  }

  double computeRadGyr(Cluster* cluster) {
    std::vector<double> com(2);
    std::set<Site*>::iterator si = cluster->members.begin();
    std::vector<double> ref = (*si)->loc;
    com[0] = 0.0;
    com[1] = 0.0;
    si++;
    for(si; si!=cluster->members.end(); si++) {
      std::vector<double> sep(2);
      sep[0] = (*si)->loc[0] - ref[0];
      sep[1] = (*si)->loc[1] - ref[1];
      mapDistance(sep);
      com[0] += sep[0];
      com[1] += sep[1];
    }
    com[0] /= cluster->members.size();
    com[1] /= cluster->members.size();
    
    si = cluster->members.begin();
    double Rgyr = 0.0;
    for(si; si!=cluster->members.end(); si++) {
      std::vector<double> sep(2);
      sep[0] = (*si)->loc[0] - ref[0];
      sep[1] = (*si)->loc[1] - ref[1];
      mapDistance(sep);
      std::vector<double> sep2(2);
      sep2[0] = sep[0] - com[0];
      sep2[1] = sep[1] - com[1];
      Rgyr += (sep2[0]*sep2[0])+(sep2[1]*sep2[1]);
    }
    Rgyr /= cluster->members.size();

    return sqrt(Rgyr);
    
  }

  void printClusterStats(std::string & clusterStatFile) {
    // order clusters by number of members before printing //

    std::multimap<int,Cluster*> orderedClusters;
    std::vector<Cluster*>::iterator cls = _clusters.begin();
    int totalSites = 0;
    for(cls; cls!=_clusters.end(); cls++) {
      orderedClusters.insert(std::pair<int,Cluster*>((*cls)->members.size(),*cls));
      totalSites += (*cls)->members.size();
    }

    std::ofstream clustFile(clusterStatFile.c_str());

    clustFile << "#total number of sites = " << totalSites << std::endl;
    clustFile << "#mass\trad\tspan\tlspan" << std::endl;

    std::multimap<int,Cluster*>::reverse_iterator ocls = orderedClusters.rbegin();
    int totalSites2 = 0;
    double corrlengthdenom = 0.0;
    double corrlengthnum = 0.0;
    double meansizenum = 0.0;
    double meansizedenom = 0.0;
    for(ocls; ocls!=orderedClusters.rend(); ocls++) {
      clustFile << ocls->first << "\t" << computeRadGyr(ocls->second);
      totalSites2 += ocls->first;
      corrlengthnum += pow(computeRadGyr(ocls->second)*(ocls->first),2.0);
      corrlengthdenom += pow(ocls->first,2.0);
      meansizenum += pow(ocls->first,2.0);
      meansizedenom += ocls->first;

      // check if cluster spans system //

      double clusSpan = 0.0;

      bool spanning = checkSpan(ocls->second,clusSpan);
      if(spanning) {
	clustFile << "\t" << 1 << "\t" << clusSpan << std::endl;
      }
      else {
	clustFile << "\t" << 0 << "\t" << clusSpan << std::endl;
      }
    }
    
    double corrlength = sqrt(2.0*corrlengthnum/corrlengthdenom);

    double meansize = meansizenum/meansizedenom;

    clustFile << std::endl << "# correlation length = " << corrlength << std::endl;
    clustFile << std::endl << "# mean cluster size = " << meansize << std::endl;

    if(totalSites != totalSites2) std::cout << "Error!" << std::endl;

    clustFile.close();

    
  }

  void findAllClusters() {
    std::vector<Site*>::iterator ni = _allsites.begin();
    for(ni; ni!=_allsites.end(); ni++) {
      if((*ni)->active && _clustersites.find(*ni) == _clustersites.end()) {
	Cluster* newCluster = new Cluster();
	addToCluster(*ni,newCluster);
	_clusters.push_back(newCluster);
      }
    }

    // check that all sites are in clusters //
  }

  void checkForSpanSite(Cluster* clus, Site* st, std::set<Site*> & visitedSites, std::vector<double> & parentSep, std::vector<double> & parentLoc, bool & span, double & largeSpan) {
    
    visitedSites.insert(st);
    
    if(span) {
      // don't do anything //
    }
    else {
      std::vector<double> sepfromParent(2);
      sepfromParent[0] = st->loc[0] - parentLoc[0];
      sepfromParent[1] = st->loc[1] - parentLoc[1];
      mapDistance(sepfromParent);
      std::vector<double> sepfromStart(2);
      sepfromStart[0] = parentSep[0] + sepfromParent[0];
      sepfromStart[1] = parentSep[1] + sepfromParent[1];
      
      // check to see if this site and the starting site span //

      double spanLength = sqrt(pow(sepfromStart[0],2)+pow(sepfromStart[1],2));

      if(spanLength > largeSpan) {
	largeSpan = spanLength;
      }
      
      if(fabs(sepfromStart[0]) >= _h[0] || fabs(sepfromStart[1]) >= _h[1]) {
	span = true;
      }
      
      // get neighbors, call function again //
      
      for(int i=0; i<st->neighbs.size(); i++) {
	Site* tmpSite = _allsites[ st->neighbs[i] ];
	if(visitedSites.find(tmpSite) == visitedSites.end() && clus->members.find(tmpSite) != clus->members.end()) {
	  checkForSpanSite(clus,tmpSite,visitedSites,sepfromStart,st->loc,span,largeSpan);
	}
      }
    }  
  }

  bool checkSpan(Cluster* clus, double & largeSpan) {
    bool span = false;
    std::set<Site*>::iterator sts = clus->members.begin();
    for(sts; sts!=clus->members.end(); sts++) {
      if(!span) {
	std::set<Site*> visitedSites;
	std::vector<double> startSep(2);
	startSep[0] = 0.0;
	startSep[1] = 0.0;
	checkForSpanSite(clus,*sts,visitedSites,startSep,(*sts)->loc,span,largeSpan);
      }
    }

    std::cout << "Analyzed span of cluster with " << clus->members.size() << " sites." << std::endl;

    return span;
  }

private:
  
  std::vector<Site*> _allsites;
  std::set<Site*> _clustersites;

  std::vector<Cluster*> _clusters;

  std::vector<double> _h;

  
};

#endif //__ClusterFinder_h__
