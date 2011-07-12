#include <iostream>
#include <iomanip>
#include <sstream>
#include <magnet/math/dilated_int.hpp>

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
  os << "\nUndilate_3 round " << i << " ( val * "
     << hexout(c::result) 
     << ") & " << hexout(z::result);
  return os.str();
}

template<size_t d, size_t i>
std::string dilate_text()
{
  typedef magnet::math::dilatedinteger::b<i,d> b;
  typedef magnet::math::dilatedinteger::y<i,d> y;

  std::ostringstream os;
  os << "\nDilate_3 round " << i << " ( val * "
     << hexout(b::result) 
     << ") & " << hexout(y::result);
  return os.str();
}

int main(int argc, char *argv[])
{
  for (size_t i(0); i < (1 << (64 / 3)); ++i)
    if (i != magnet::math::undilate<3>(magnet::math::dilate<3>(i)))
      std::cout << "\nError for a value of i=" << i << " " << magnet::math::dilate<3>(i);

  for (size_t i(0); i < (1 << (64 / 4)); ++i)
    if (i != magnet::math::undilate<4>(magnet::math::dilate<4>(i)))
      return -1;
      
  for (size_t i(0); i < (1 << (64 / 5)); ++i)
    if (i != magnet::math::undilate<5>(magnet::math::dilate<5>(i)))
      return -1;

  for (size_t i(0); i < (1 << (64 / 6)); ++i)
    if (i != magnet::math::undilate<6>(magnet::math::dilate<6>(i)))
      return -1;

  for (size_t i(0); i < (1 << (64 / 7)); ++i)
    if (i != magnet::math::undilate<7>(magnet::math::dilate<7>(i)))
      return -1;

  for (size_t i(0); i < (1 << (64 / 8)); ++i)
    if (i != magnet::math::undilate<8>(magnet::math::dilate<8>(i)))
      return -1;

  for (size_t i(0); i < (1 << (64 / 9)); ++i)
    if (i != magnet::math::undilate<9>(magnet::math::dilate<9>(i)))
      return -1;

  for (size_t i(0); i < (1 << (64 / 10)); ++i)
    if (i != magnet::math::undilate<10>(magnet::math::dilate<10>(i)))
      return -1;

  return 0;
}
