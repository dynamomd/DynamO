#include "dilatedint.hpp"

const unsigned int MaskedInteger::Smask = (((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S);
const unsigned int MaskedInteger::maxVal = (((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S);
const unsigned int MaskedInteger::dilatedMaxVal = ((((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S * 3)) & mask; 

