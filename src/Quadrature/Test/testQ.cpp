#include "LineQuadrature.h"
#include <iostream>

using namespace voom;

int main()
{
  LineQuadrature q1(1);
  q1.check(1);

  LineQuadrature q3(3);
  q3.check(3);

  LineQuadrature q5(5);
  q5.check(5);

  return 0;
}



