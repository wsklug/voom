#include "QuadQuadrature.h"
#include <iostream>

using namespace std;
using namespace voom;

int main()
{
  for(int i=1; i<=3; i++) {
    QuadQuadrature quad(i);
    
    quad.check(i);
  }  

  return 0;
}
