#include "Hermite.h"
#include <iostream>
//#include <rand>

using namespace std;
using namespace voom;

int main()
{
	// check 
  
	srand(time(0));
	double xi = 2.0*(double)(rand())/RAND_MAX-1.0;
	Hermite shp(xi);

	int n = shp.derivatives().size();
	std::cout << n << " nodes." << std::endl;

	Hermite::DerivativeContainer DN;
	DN.resize(shp.derivatives().size());

	Hermite::DerivativeContainer DDN;
	DDN.resize(shp.secondDerivatives().size());

	double dxi = 1.0e-8;
	
	shp.compute(xi+dxi);
	for(int i=0; i<4; i++) {
	  DN[i] = shp.functions()[i];
	  DDN[i] = shp.derivatives()[i];
	}

	shp.compute(xi-dxi);
	for(int i=0; i<4; i++) {
	  DN[i] -= shp.functions()[i];
	  DDN[i] -= shp.derivatives()[i];
	}

	shp.compute(xi);
	for(int i=0; i<4; i++) {
	  DN[i] /= 2.0*dxi;
	  DDN[i] /= 2.0*dxi;
	}

	double error = 0.0;
	double norm = 0.0;
	double errorD = 0.0;
	double normD = 0.0;
	for(int a=0; a<shp.derivatives().size(); a++) {
	  double e=DN[a]-shp.derivatives()[a];
	  error += e*e;
	  norm += shp.derivatives()[a]*shp.derivatives()[a];
	  e=DDN[a]-shp.secondDerivatives()[a];
	  errorD += e*e;
	  normD += shp.secondDerivatives()[a]*shp.secondDerivatives()[a];
	}
	norm = sqrt(norm);
	error = sqrt(error);
	
	normD = sqrt(normD);
	errorD = sqrt(errorD);
	
	std::cout << "Error = " << error << " Norm = " << norm << std::endl;
	std::cout << "ErrorD = " << errorD << " NormD = " << normD << std::endl;

	double tol = 1.0e-6;

	if( std::abs(error) < tol*norm && std::abs(errorD) < tol*normD  ) {
	  std::cout<<"Shape consistency check PASSED!"<<std::endl;
	  return 0;
	}

	std::cout<<"Shape consistency check FAILED!"<<std::endl
		 << std::setw(10) << "n"
		 << std::setw(24) << "analytical"
		 << std::setw(24) << "numerical" << std::endl;
	for(int a=0; a<shp.derivatives().size(); a++) {
	    std::cout << std::setw(10) << a 
		      << std::setw(24) << shp.derivatives()[a]
		      << std::setw(24) << DN[a]
		      << std::endl;	
	}

	return 0;
}
