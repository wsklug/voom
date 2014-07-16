#include "ClusterFinder.h"

#include<string>
#include<set>
#include<vector>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char* argv[]) {
  char* connFileName = argv[1];
  std::string connFN(connFileName);

  std::vector<double> ssize(2);
  ssize[0] = atof(argv[2]);
  ssize[1] = atof(argv[3]);

  ClusterFinder* cfinder = new ClusterFinder(connFN,ssize);

  std::cout << "Successfully read in cluster information." << std::endl;

  if(argc > 4) {
    std::string activityFile(argv[4]);
    cfinder->readActivity(activityFile);

    std::cout << "Finished reading in state of segments from " << activityFile << std::endl;

  }

  cfinder->findAllClusters();

  std::cout << "Please enter the name of a file in which to store cluster data: ";

  char clustFile[256];
  
  std::cin >> clustFile;
  
  std::string clustFN(clustFile);

  cfinder->printClusterStats(clustFN);


}
