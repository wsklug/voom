#include "LoopShellShape.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"
#include "ShapeQ4.h"
#include "ShapeTet10.h"
#include <iostream>
//#include <rand>

using namespace std;
using namespace voom;

int main()
{
  srand(time(0));
  
  ShapeTet10::CoordinateArray s;
  s(0) = 1.0*rand()/RAND_MAX;
  s(1) = (1.0-s(0))*rand()/RAND_MAX;
  s(2) = (1.0-s(0))*rand()/RAND_MAX;
  
  ShapeTet10 tet10(s);
  tet10.checkConsistency(s);
  
  return 0;
}
