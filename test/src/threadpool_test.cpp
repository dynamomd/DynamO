#include <iostream>
#include <vector>
#include <stdexcept>
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

void function2(int i)
{  
  float sum = 0;
  for (int j = 0; j < i; ++j)
    sum *= j;

  sums[i] = sum;
}

struct A
{
  void memberFunc() { std::cerr << "Inside memberfunc\n"; }

  void memberFunc2(int i) { std::cerr << "Inside memberfunc2, i=" << i << "\n"; }

  void memberFunc3(int& i, int j) 
  { std::cerr << "Inside memberfunc3, i=" << i << ", j=" << j << "\n"; }
};

int main()
{
  int N = 1000;
  sums.resize(N);

  A Aclass;

  magnet::thread::ThreadPool pool;

  pool.setThreadCount(4);

  std::cerr << "Using " << pool.getThreadCount() << " threads\n";
  
  //const int& val1 = 2;
  //const int& val2 = 2;
  int val3 = 2;
  int& val4 = val3;
  
  pool.queueTask(magnet::function::Task::makeTask(&A::memberFunc, &Aclass));
  pool.queueTask(magnet::function::Task::makeTask(&A::memberFunc2, &Aclass, 2));
  pool.queueTask(magnet::function::Task::makeTask(&A::memberFunc3, &Aclass, val4, 4));


  for (size_t loop(0); loop < 1000; ++loop)
    {
      //std::cerr << "Function 1\n";
      
      for (int i = 0; i < N; ++i)   
	pool.queueTask(magnet::function::Task::makeTask(function1, i));
      
      pool.wait();
      
      for (int i = N-1; i >= 0; --i)
	{
	  float tmp = sums[i];
	  function1(i);
	  if (sums[i] != tmp) 
	    {
	      std::cerr << "Failure in loop " << loop << "\n";
	      throw std::runtime_error("Muck up in function 1");
	    }
	}
      
      //std::cerr << "Function 2\n";
      for (int i = 0; i < N; ++i)    
	pool.queueTask(magnet::function::Task::makeTask(function2, i));
      
      pool.wait();
      
      for (int i = N-1; i >= 0; --i)
	{
	  float tmp = sums[i];
	  function2(i);
	  if (sums[i] != tmp)
	    {
	      std::cerr << "Failure in loop " << loop << "\n";
	      throw std::runtime_error("Muck up in function 2");
	    }
	}
    }

  std::cerr << "Finished\n";

  return 0;
}
