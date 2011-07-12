#include <magnet/math/dilated_int.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

std::string hexout(size_t val)
{
  std::ostringstream os;
  os << "0x"<< std::hex << std::setw(sizeof(size_t)*2) << std::setfill ('0') << val;
  return os.str();
}

template<size_t d, size_t i>
std::string undilate_text()
{
  typedef magnet::math::dilatedinteger::c<i,d> c;
  typedef magnet::math::dilatedinteger::z<i,d> z;

  std::ostringstream os;
  os << "\nUndilate<" << d <<"> round " << i << " ( val * " << hexout(c::result) << ") & " << hexout(z::result);
  return os.str();
}

template<size_t d, size_t i>
std::string dilate_text()
{
  typedef magnet::math::dilatedinteger::b<i,d> b;
  typedef magnet::math::dilatedinteger::y<i,d> y;
 
  std::ostringstream os;
  os << "Dilate<" << d << "> round " << i << " ( val * " << hexout(b::result) << ") & " << hexout(y::result)
     << "\n";
  return os.str();
}

template<size_t i>
std::string dilate_text_2()
{
  static size_t shiftval = magnet::math::dilatedinteger::dilate<2>::shiftval<i>::result;

  std::ostringstream os;
  os << "Dilate<2> round " << i << " ( val | (val <<  " << shiftval << ")) & " 
     << hexout(magnet::math::dilatedinteger::y<i,2>::result)
     << "\n";
  return os.str();
}

template<size_t d>
void test()
{
  size_t maxval = magnet::math::dilatedinteger::maxDilatableValue<d>::result;
  std::cout << "Testing d = " << d << " dilation, max val = " << maxval << std::endl;
  
  size_t testvals = maxval;
  if (testvals > 2097151) {testvals = 2097151; std::cout << "   limiting to 2097151 tests\n"; }
  for (size_t i(0); i < testvals; ++i)

    if (i != magnet::math::undilate<d>(magnet::math::dilate<d>(i)))
      {
	std::cerr << "Failure for i = " << i << "\n";
	throw std::runtime_error("Undilate and dilate are not symmetric");
      }

  //Double check the max value, just to make sure nothing edge-casey
  //is going on
  if (maxval != magnet::math::undilate<d>(magnet::math::dilate<d>(maxval)))
    {
      std::cerr << "Failure for i = maxval\n";
      throw std::runtime_error("Undilate and dilate are not symmetric for the max value");
    }
}

int main(int argc, char *argv[])
{
  try {
    test<2>();
    test<3>();
    test<4>();
    test<5>();
    test<6>();
    test<7>();
    test<8>();
    test<9>();
    test<10>();
  } catch (std::exception&)
    {
      std::cerr << "Failed the dilation test";
    }
}
