#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <unistd.h>
#include <math.h>

using namespace std;

/*! This program reads in the MTiteration.vtk files
and return number of hexamers with direction +pi/3, 0.0 and -pi/3
 */
int main(int argc, char* argv[])
{
  // Input parameters
  unsigned int Nmax = 0;  // Last MT iteration
  unsigned int Nhex = 0;  // Number of hexamer
  unsigned int Nt = 0;    // Number of triangles per hexamer

  // Open input parameters file and read parameters
  if(argc < 2){ 
    cout << "Input file missing." << endl; 
    return(0);
  }
  string parameterFileName = argv[1];
  ifstream inp;
 
  inp.open(parameterFileName.c_str(), ios::in);
  if (!inp) {
    cout << "Cannot open input file: " << parameterFileName << endl;
    return(0);
  }
  inp >> Nmax >> Nhex >> Nt;
  inp.close();

  // List input parameters
  cout << " Last MT iteration               : " << Nmax << endl
       << " Number of hexamer               : " << Nhex << endl
       << " Number of triangles per hexamer : " << Nt   << endl;
  
  // Number of elements in the hexamers
  unsigned int Nel = Nt*Nhex;



  // Loop over MT output files - write results in SatesSummary.dat
  ofstream StatesSummary("StatesSummary.dat");
  if (!StatesSummary) {
    std::cout << "Cannot open output StatesSummary.dat" << std::endl;
    return(0);
  }	
  StatesSummary  << "MTiter Nred Nwhite Nblue" << endl;
	 


  string token;
  double angle = 0.0, alpha = M_PI/3.0, tol = 1.0e-4;
  unsigned int iter = 0, j = 0;
  unsigned int Nred = 0, Nwhite = 0, Nblue = 0;
  ifstream ifs;
  for (iter = 0; iter < Nmax+1; iter++)
  {
    Nred = 0;
    Nwhite = 0;
    Nblue = 0;
    stringstream MToutput;
    MToutput << "MTiteration" << iter << ".vtk";
    ifs.open( (MToutput.str()).c_str(), ios::in);
    if (!ifs) {
      cout << "Cannot open input file: " << MToutput << endl;
      return(0);
    }

    while( token != "StretchDirection" ) 
    {
      ifs >> token;
    }
    ifs >> token;
    ifs >> token;
    ifs >> token;
    ifs >> token;
    for (j = 0; j < Nel; j++)
    { 
      ifs >> angle;
      if (fabs(angle - alpha) < tol)
      {
	Nred++;
      }
      else if (fabs(angle + alpha) < tol)
      {
	Nblue++;
      } 
      else if (fabs(angle) < tol)
      {
	Nwhite++;
      } 
      else
      {
	cout << "Error ! Angle = " << angle << endl;
      }
    }
    ifs.close();

    if ( (Nred + Nwhite + Nblue)/Nt != Nhex)
    {
      cout << "Error ! N hexamers = " << (Nred + Nwhite + Nblue)/Nt << endl;
    }
    StatesSummary  << iter << " " << Nred/Nt << " " << Nwhite/Nt << " " << Nblue/Nt << endl;
  }

  StatesSummary.close();

  return (0);
  
}


