#include "HelperFunctions.h"

using namespace tvmet;
using namespace std;
using namespace voom;

int main(int argc, char* argv[]) {
	std::stringstream sstm;
	std::string rName;
	sstm.str("");
	sstm.clear();
	std::vector<std::string> fileNames;
	for (int i = 0; i < 547; i++) {
		sstm << "T7-relaxed-" << i << ".vtk";
		rName = sstm.str();
		fileNames.push_back(rName);
		sstm.str("");
		sstm.clear();
	}
	writeEdgeStrainVtk(fileNames, 1.001137586, 10.0);
	return 0;
}
