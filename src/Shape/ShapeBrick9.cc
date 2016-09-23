#include "ShapeBrick9.h"

namespace voom {  
  void ShapeBrick9::compute(const CoordinateArray & s) {

    // Consistent with node numbering in ShapeBrick9.h
    _functions[0] = (1.0-s(0)-s(1))*(1.0+s(2)*(2*s(2)-3.0));
    _functions[1] = s(0)*(1.0+s(2)*(2*s(2)-3.0));
    _functions[2] = s(1)*(1.0+s(2)*(2*s(2)-3.0));

    _functions[3] = (1.0-s(0)-s(1))*4.0*s(2)*(1.0-s(2));
    _functions[4] = s(0)*4.0*s(2)*(1.0-s(2));
    _functions[5] = s(1)*4.0*s(2)*(1.0-s(2));

    _functions[6] = (1.0-s(0)-s(1))*s(2)*(2.0*s(2)-1.0);
    _functions[7] = s(0)*s(2)*(2.0*s(2)-1.0);
    _functions[8] = s(1)*s(2)*(2.0*s(2)-1.0);

    _derivatives[0] = s(2)*(3.0-2*s(2))-1.0, s(2)*(3.0-2*s(2))-1.0, (1.0-s(0)-s(1))*(4.0*s(2)-3.0);
    _derivatives[1] = (1.0+s(2)*(2*s(2)-3.0)), 0.0, s(0)*(4.0*s(2)-3.0);
    _derivatives[2] = 0.0, (1.0+s(2)*(2*s(2)-3.0)), s(1)*(4.0*s(2)-3.0);

    _derivatives[3] = -4.0*s(2)*(1.0-s(2)), -4.0*s(2)*(1.0-s(2)), (1.0-s(0)-s(1))*(4.0-8.0*s(2));
    _derivatives[4] = 4.0*s(2)*(1.0-s(2)), 0.0, s(0)*(4.0-8.0*s(2));
    _derivatives[5] = 0.0, 4.0*s(2)*(1.0-s(2)), s(1)*(4.0-8.0*s(2));

    _derivatives[6] = s(2)*(1.0-2.0*s(2)), s(2)*(1.0-2.0*s(2)), (1.0-s(0)-s(1))*(4.0*s(2)-1.0);
    _derivatives[7] = s(2)*(2.0*s(2)-1.0), 0.0, s(0)*(4.0*s(2)-1.0);
    _derivatives[8] = 0.0, s(2)*(2.0*s(2)-1.0), s(1)*(4.0*s(2)-1.0);
    
  }
};

