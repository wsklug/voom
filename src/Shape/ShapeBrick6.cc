#include "ShapeBrick6.h"

namespace voom {  
  void ShapeBrick6::compute(const CoordinateArray & s) {

    // Consistent with node numbering in ShapeBrick6.h
    _functions[0] = (1.0-s(0)-s(1))*(1.0-s(2));
    _functions[1] = s(0)*(1.0-s(2));
    _functions[2] = s(1)*(1.0-s(2));
    _functions[3] = (1.0-s(0)-s(1))*s(2);
    _functions[4] = s(0)*s(2);
    _functions[5] = s(1)*s(2);

    _derivatives[0] = s(2)-1.0, s(2)-1.0, s(0)+s(1)-1.0;
    _derivatives[1] = 1.0-s(2), 0.0, -s(0);
    _derivatives[2] = 0.0, 1.0-s(2), -s(1);
    _derivatives[3] = -s(2), -s(2), 1.0-s(0)-s(1);
    _derivatives[4] = s(2), 0.0, s(0);
    _derivatives[5] = 0.0, s(2), s(1);
    
  }
};

