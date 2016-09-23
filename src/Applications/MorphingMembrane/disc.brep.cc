#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

int main() {
  int npts=24;

  char fname[50];
  sprintf(fname,"disc-%d.brep",npts);
  std::ofstream ofs(fname);
  ofs.setf(ios::showpoint);
  ofs.setf(ios::fixed);
//   ofs.precision(8);
  ofs << npts << std::endl;

  double r[npts][2];
  for(int i=0; i<npts; i++) {
    double theta = 2.0*M_PI*i/npts;
    r[i][0] = cos(theta);
    r[i][1] = sin(theta);
    ofs << i+1 
	<< '\t' << r[i][0]
	<< '\t' << r[i][1] << std::endl;
  }

  // print edges
  ofs << npts << std::endl;
  for(int i=0; i<npts; i++) {
    int iA=i;
    int iB=(i+1<npts) ? i+1 : 0;
    ofs << i+1 << "\t"
	<< iA+1 << "\t" << iB+1
	<< "\t" << r[iA][0]
	<< "\t" << r[iA][1]
	<< "\t" << r[iB][0]
	<< "\t" << r[iB][1]
	<< std::endl;
  }

  // print loops
  ofs << 1 << std::endl
      << "1 " << npts << std::endl;
  for(int i=0; i<npts; i++) {
    ofs << i+1 << " ";
  }
  ofs << std::endl;

  // print subbody and body
  ofs << 1 << std::endl
      << 1 << " " << 1 << std::endl
      << 1 << std::endl
      << 1 << " " << 1 << std::endl
      << 1 << std::endl;
  
  return 0;
}
