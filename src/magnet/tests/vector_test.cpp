#include <magnet/math/vector.hpp>
#include <iostream>
#include <cmath>

bool err(double val, double expected)
{
  return std::abs(val / expected - 1) > 0.0001;
}

int main()
{
  Vector A(0,0,0);
  Vector B(1,1,1);

  if (err((A+B).nrm(), std::sqrt(3)))
    { std::cout << "(A+B).nrm() is wrong"; return 1; }

  return 0;
}
