#include <iostream>
#include <vector>
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

struct A
{
  void memberFunc() { std::cerr << "Inside memberfunc\n"; }

  void memberFunc2(int i) { std::cerr << "Inside memberfunc2, i=" << i << "\n"; }

  void memberFunc3(int& i, int j) 
  { std::cerr << "Inside memberfunc2, i=" << i << ", j=" << j << "\n"; }
};

int main()
{
  size_t N = 8000;
  sums.resize(N);

  A Aclass;

  ThreadPool pool;

  pool.setThreadCount(3);
  
  const int& val1 = 2;
  const int& val2 = 2;
  int val3 = 2;
  int& val4 = val3;
  

  pool.invoke(ThreadPool::makeTask(&A::memberFunc, &Aclass));
  
  pool.invoke(ThreadPool::makeTask(&A::memberFunc2, &Aclass, 2));
  
  pool.invoke(ThreadPool::makeTask<void, A, int&, int>(&A::memberFunc3, &Aclass, val4, 4));

  for (int i = 0; i < N; ++i)    
    pool.invoke(ThreadPool::makeTask(function1, i+0));

  std::cerr << "Entering Wait\n";
  pool.wait(); 
  std::cerr << "Finished\n";

  return 1;
}
