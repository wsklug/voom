#include "../LoopShellShape.h"
#include "../ShapeTri3.h"
#include "../ShapeTri6.h"
#include "../ShapeQ4.h"
#include "../ShapeBrick6.h"
#include "../ShapeBrick9.h"
#include <iostream>
//#include <rand>

using namespace std;
using namespace voom;

int main()
{
  srand(time(0));
  
  ShapeBrick6::CoordinateArray s;
  s(0) = 1.0*rand()/RAND_MAX;
  s(1) = (1.0-s(0))*rand()/RAND_MAX;
  s(2) = (1.0-s(0))*rand()/RAND_MAX;

  ShapeBrick6 Brick6(s);
  Brick6.checkConsistency(s);
  Brick6.checkC0Completeness();



  ShapeBrick9 Brick9(s);
  Brick9.checkConsistency(s);
  Brick9.checkC0Completeness();
  
  return 0;
}
