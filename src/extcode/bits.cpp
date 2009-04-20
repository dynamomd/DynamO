#include <iostream>
#include <limits>
#include <assert.h>

//bit masks for the dimensions
const unsigned int dimbits[] = {0x49249249, 0x92492492, 0x24924924};

inline unsigned int incrementMorton(const unsigned int& MortonNum, const unsigned int& mask)
{
  return (MortonNum & (dimbits[mask] ^ 0xFFFFFFFF))
    | ((MortonNum & dimbits[mask]) - dimbits[mask]) & dimbits[mask];
}

inline unsigned int decrementMorton(const unsigned int& MortonNum, const unsigned int& mask)
{
  return (MortonNum & (dimbits[mask] ^ 0xFFFFFFFF))
    | ((MortonNum & dimbits[mask]) - 1) & dimbits[mask];
}

/*inline unsigned int addMorton(const unsigned int& MortonNum, const unsigned int& add,
			      const unsigned int& mask)
{
  return (MortonNum & (dimbits[mask] ^ 0xFFFFFFFF))
    | ((MortonNum & dimbits[mask]) - 1) & dimbits[mask];
    }*/

unsigned int get3DMortonNum(const unsigned char& x,
			    const unsigned char& y,
			    const unsigned char& z)
{
  assert(sizeof(unsigned int) == 4);
  unsigned int coord(0);
  
  for (int i = 0; i < sizeof(x) * CHAR_BIT; i++)
    {
      coord |=
	((x & (1 << i)) << (i*2))
	| ((y & (1 << i)) << (i*2 + 1))
	| ((z & (1 << i)) << (i*2 + 2));
    }

  return coord;
}

void get3DMortonCoords(const unsigned int& MortonNum)
{
  assert(sizeof(unsigned int) == 4);
 
  unsigned int x(0), y(0), z(0);

  for (int i = 0; i < sizeof(unsigned char) * CHAR_BIT; i++)
    {
      x |= (MortonNum & (1 << i*3)) >> (i*2);
      y |= (MortonNum & (1 << i*3+1)) >> (i*2+1);
      z |= (MortonNum & (1 << i*3+2)) >> (i*2+2);
    }

  std::cout //<< std::hex << std::showbase 
	    << x << " " << y << " " << z << std::endl;
}

int main()
{
  unsigned int MortonNum = get3DMortonNum(2,0,3);
  get3DMortonCoords(MortonNum);
  get3DMortonCoords(decrementMorton(MortonNum,0));
  get3DMortonCoords(decrementMorton(MortonNum,1));
  get3DMortonCoords(decrementMorton(MortonNum,2));
  
}
