#include "BrickQuadrature.h"
#include <iostream>

using namespace std;
using namespace voom;

int main()
{
  BrickQuadrature BrickQuad1(1.0);  
  for (unsigned int i=0; i<=1; i++)
  {
    for (unsigned int j=0; j<=1; j++)
    {
      BrickQuad1.check(i, j);
      std::cout << endl;
    }
  }   

  BrickQuadrature BrickQuad2(2.0);  
  for (unsigned int i=0; i<=2; i++)
  {
    for (unsigned int j=0; j<=1; j++)
    {
      BrickQuad2.check(i, j);
      std::cout << endl;
    }
  }   
  
  BrickQuadrature BrickQuad3(3.0);  
  for (unsigned int i=0; i<=2; i++)
  {
    for (unsigned int j=0; j<=2; j++)
    {
      BrickQuad3.check(i, j);
      std::cout << endl;
    }
  }   
  
  

  return 0;
}
