#include <iostream>
#include <boost/bind.hpp>
#include "../src/extcode/threadpool.hpp"

std::vector<float> sums;

void function1(int i) 
{  
  float sum = 0;
  for (int j = 0; j < i; ++j)
    sum += j;

  sums[i] = sum;
}

int main()
{
  size_t N = 40000;
  sums.resize(N);

  ThreadPool pool;

  pool.setThreadCount(3);

  for (int i = 0; i < N; ++i)
    pool.invoke(boost::bind(&function1, i));

  std::cerr << "Entering Wait\n";
  pool.wait(); 
  std::cerr << "Finished\n";

  return 1;
}
