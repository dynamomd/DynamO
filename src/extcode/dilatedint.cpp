#include "dilatedint.hpp"

const unsigned int MaskedInteger::Smask = (((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S);
const unsigned int MaskedInteger::maxVal = (((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S);
const unsigned int MaskedInteger::dilatedMaxVal = ((((unsigned int)0) - 1) >> (sizeof(unsigned int)*CHAR_BIT - S * 3)) & mask; 

/*

int main()
{
  dilatedCoords test(1023,3,0);
  
  std::cout << test.getMortonNum() << std::endl;
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;

  ++test.data[0];
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;

  --test.data[1];
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal()
	    << " " << test.data[2].getRealVal() <<  std::endl;

  test.data[1] = test.data[1] + MI(5);
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;

  --test.data[2];
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;

  test = dilatedCoords(153391707);
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;

  test.data[0] = 21;
  std::cout << test.data[0].getRealVal() << " " << test.data[1].getRealVal() 
	    << " " << test.data[2].getRealVal() <<  std::endl;
  
  std::cout << std::endl;

  print_bits(MI(1).getDilatedVal(),std::cout);
  std::cout << std::endl;

  print_bits((MI(2)).getDilatedVal(),std::cout);
  std::cout << std::endl;

  MI i(1);
  print_bits((i+2).getDilatedVal(), std::cout);
  std::cout << std::endl;

  std::cout << MI::maxVal << std::endl;
  print_bits(MI::dilatedMaxVal, std::cout);
  std::cout << std::endl;
}
*/
