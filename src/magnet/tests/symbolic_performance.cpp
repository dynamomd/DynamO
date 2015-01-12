/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <chrono>
#include <cmath>
#include <random>
#include <magnet/math/symbolic.hpp>

std::mt19937 RNG;
auto dist = std::uniform_real_distribution<double>(-1, 1);

struct TimeScope {
  TimeScope(std::string message): _message(message), _counter(0) {
    _start = std::chrono::high_resolution_clock::now();
  }
  
  void operator++() {_counter++;}
  
  ~TimeScope() {
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration <double, std::nano>(stop - _start).count() / double(_counter);
    std::cout << _message << " " << duration << " ns/call" << std::endl;
  }

  std::string _message;
  std::chrono::high_resolution_clock::time_point _start;
  size_t _counter;
};

template<class T>
bool check_close(T f1, T f2) {
  return std::abs(std::abs(f1) - std::abs(f2)) > (std::abs(f1) + std::abs(f2)) * 1e-8;
}

double std_val, psym_val, sym_val;

void testValues() {
  if (check_close(sym_val, std_val) || check_close(psym_val, std_val)) {
    std::cout << "  WARNING! Mismatch in results!" << std::endl;
    std::cout << "   standard                = " << std_val << std::endl;
    std::cout << "   pre-calculated symbolic = " << psym_val << std::endl;	      
    std::cout << "   symbolic                = " << sym_val << std::endl;	      
  }
}

int main(int argv, const char** argc)
{
  const size_t tests = 1000000;
  //Just for convenience, define x
  magnet::math::Variable<'x'> x;

  /////////////////////////////////////////////////////////////////
  //////////////////////////// TEST ///////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "±" << std::endl;
  std::cout << "f(x) = " << expand(x*x+2*x-3) << std::endl;  
  RNG.seed(12345);
  {
    TimeScope timer(" Standard                 ");
    std_val = 0;
    for (size_t i(0); i < tests; ++i) {
      auto y = dist(RNG);
      std_val += y*y + 2*y - 3;
      ++timer;
    }
  }

  RNG.seed(12345);
  {
    TimeScope timer(" Symbolic (pre-calculated)");
    using namespace magnet::math;
    auto f = x*x+2*x-3;
    psym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      psym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }

  RNG.seed(12345);
  {
    TimeScope timer(" Symbolic (pc & expanded) ");
    using namespace magnet::math;
    auto f = expand(x*x+2*x-3);
    psym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      psym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }

  RNG.seed(12345);
  {
    TimeScope timer(" Symbolic                 ");
    using namespace magnet::math;
    sym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      sym_val += substitution(x*x+2*x-3, x==dist(RNG));
      ++timer;
    }
  }

  testValues();

  /////////////////////////////////////////////////////////////////
  //////////////////////////// TEST ///////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "\nf(x) = sin(x^2 + 2 x - 3) - 2 * cos(x)" << std::endl;  
  RNG.seed(12345);
  {
    TimeScope timer(" Standard                 ");
    std_val = 0;
    for (size_t i(0); i < tests; ++i) {
      auto y = dist(RNG);
      std_val += std::sin(y*y + 2*y - 3) - 2*std::cos(y);
      ++timer;
    }
  }

  RNG.seed(12345);
  {
    TimeScope timer(" Symbolic (pre-calculated)");
    using namespace magnet::math;
    psym_val = 0;
    auto f = sin(x*x+2*x-3) - 2*cos(x);
    for (size_t i(0); i < tests; ++i) {
      psym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }

  RNG.seed(12345);
  {
    TimeScope timer(" Symbolic                 ");
    using namespace magnet::math;
    sym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      sym_val += substitution(sin(x*x+2*x-3) - 2*cos(x), x==dist(RNG));
      ++timer;
    }
  }
  testValues();

  if (check_close(sym_val, std_val)) {
    std::cout << "  WARNING! Mismatch in results!" << std::endl;
    std::cout << "   standard = " << std_val << std::endl;
    std::cout << "   symbolic = " << sym_val << std::endl;	      
  }

  /////////////////////////////////////////////////////////////////
  //////////////////////////// TEST ///////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "\nf'(x), where f(x) = sin(x^2 + 2 x - 3) - 2 * cos(x)" << std::endl;  
  {
    RNG.seed(12345);
    TimeScope timer(" Standard                 ");
    std_val = 0;
    for (size_t i(0); i < tests; ++i) {
      auto y = dist(RNG);
      std_val += (2 * y + 2) * std::cos(y*y+2*y-3) + 2 * std::sin(y);
      ++timer;
    }
  }

  {
    RNG.seed(12345);
    TimeScope timer(" Symbolic (pre-calculated)");
    psym_val = 0;
    auto f = derivative(sin(x*x+2*x-3) - 2*cos(x), x);
    using namespace magnet::math;
    for (size_t i(0); i < tests; ++i) {
      psym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }

  {
    RNG.seed(12345);
    TimeScope timer(" Symbolic                 ");
    sym_val = 0;
    auto f = derivative(sin(x*x+2*x-3) - 2*cos(x), x);
    using namespace magnet::math;
    for (size_t i(0); i < tests; ++i) {
      sym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }
  testValues();

  /////////////////////////////////////////////////////////////////
  //////////////////////////// TEST ///////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "\n5th order Taylor expansion of f(x) = sin(x^2 + 2 x - 3) - 2 * cos(x)" << std::endl;  
  {
    RNG.seed(12345);
    TimeScope timer(" Standard                 ");
    std_val = 0;
    for (size_t i(0); i < tests; ++i) {
      auto y = dist(RNG);
      auto y2 = y*y;
      auto y3 = y2*y;
      auto y4 = y3*y;
      auto y5 = y4*y;
      std_val += 181.2677681603864*y5-2837.24417459026*y4+17665.06763284448*y3-54699.69647619253*y2+84257.31059283158*y-51661.33568865078;
      ++timer;
    }
  }

  {
    using namespace magnet::math;
    auto f = taylor_series<5, 'x'>(sin(x*x + 2*x - 3) - 2 * cos(x), 3.0);
    RNG.seed(12345);
    TimeScope timer(" Symbolic (pre-calculated)");
    psym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      psym_val += substitution(f, x==dist(RNG));
      ++timer;
    }
  }

  {
    using namespace magnet::math;
    RNG.seed(12345);
    TimeScope timer(" Symbolic                 ");
    sym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      sym_val += substitution(taylor_series<5, 'x'>(sin(x*x + 2*x - 3) - 2 * cos(x), 3.0), x==dist(RNG));
      ++timer;
    }
  }

  testValues();

  /////////////////////////////////////////////////////////////////
  //////////////////////////// TEST ///////////////////////////////
  /////////////////////////////////////////////////////////////////
  std::cout << "\nSolve roots of f(x) = x^2 + 2 x - 3" << std::endl;  
  {
    RNG.seed(12345);
    TimeScope timer(" Standard                 ");
    std_val = 0;
    for (size_t i(0); i < tests; ++i) {
      double root1 = dist(RNG);
      double root2 = dist(RNG);
      double a = dist(RNG);
      double b = a * (-root1 -root2);
      double c = a * root1 * root2;
      double arg = b*b-4*a*c;
      double root1_solved = (-b +std::sqrt(arg)) / (2 * a);
      double root2_solved = (-b -std::sqrt(arg)) / (2 * a);
      std_val += root1_solved + root2_solved;
      ++timer;
    }
  }

  {
    using namespace magnet::math;
    RNG.seed(12345);
    TimeScope timer(" Symbolic (Polynomial)");
    psym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      const double root1 = dist(RNG);
      const double root2 = dist(RNG);
      const double a = dist(RNG);
      const double b = a * (-root1 -root2);
      const double c = a * root1 * root2;
      auto f = Polynomial<2>{c,b,a};
      auto roots = solve_real_roots(f);
      psym_val += roots[0] + roots[1];
      ++timer;
    }
  }

  {
    using namespace magnet::math;
    RNG.seed(12345);
    TimeScope timer(" Symbolic (expanded)      ");
    sym_val = 0;
    for (size_t i(0); i < tests; ++i) {
      const double root1 = dist(RNG);
      const double root2 = dist(RNG);
      const double a = dist(RNG);
      auto f = expand(a * (x - root1) * (x - root2));
      auto roots = solve_real_roots(f);
      sym_val += roots[0] + roots[1];
      ++timer;
    }
  }

  testValues();

}
