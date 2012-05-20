/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#include <magnet/CL/detail/common.hpp>

#define STRINGIFY(A) #A

namespace magnet {
  namespace CL {
    template<class T>
    class heapSort : public detail::Program
    {
      cl::Kernel _sortKernel, _dataSortKernel;
    
    public:
      void build(cl::CommandQueue queue, cl::Context context)
      {
	Program::build(queue, context);

	// set up kernel
	_sortKernel = cl::Kernel(_program, "heapSort");
	_dataSortKernel = cl::Kernel(_program, "heapSortData");	
      }

      void operator()(cl::Buffer input, cl_uint ascending = true)
      {
	const cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);

	cl::KernelFunctor clsort = _sortKernel.bind(_queue, cl::NDRange(1), 
						     cl::NDRange(1));      
	clsort(input, size);
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput)
      {
	const cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	if (size != dataInput.getInfo<CL_MEM_SIZE>() / sizeof(cl_uint))
	  M_throw() << "Data-key buffer size mismatch";
	

	cl::KernelFunctor clsort = _dataSortKernel.bind(_queue, cl::NDRange(1), 
							cl::NDRange(1));      
	clsort(keyInput, dataInput, size);

      }

      virtual std::string initKernelSrc()
      {
	return std::string("#define keyType ") + detail::traits<T>::kernel_type() + "\n"
	  STRINGIFY(
void siftDown(__global keyType* numbers, int root, int bottom)
{
  int done = 0;
 
  done = 0;
  while ((root*2 <= bottom) && (!done))
    {
      int maxChild = root * 2 + 1;

      if (root*2 == bottom)
	maxChild = root * 2;
      else if (numbers[root * 2] > numbers[root * 2 + 1])
	maxChild = root * 2;
    
      if (numbers[root] < numbers[maxChild])
	{
	  keyType temp = numbers[root];
	  numbers[root] = numbers[maxChild];
	  numbers[maxChild] = temp;
	  root = maxChild;
	}
      else
	done = 1;
    }
}

__kernel void heapSort(__global keyType* numbers, uint array_size)
{
  int i;
 
  for (i = (array_size / 2) - 1; i >= 0; i--)
    siftDown(numbers, i, array_size - 1);
 
  for (i = array_size-1; i >= 1; i--)
    {
      keyType temp = numbers[0];
      numbers[0] = numbers[i];
      numbers[i] = temp;
      siftDown(numbers, 0, i-1);
    }
}

void siftDownData(__global keyType* numbers, __global uint* data,  int root, int bottom)
{
  int done = 0;
 
  done = 0;
  while ((root*2 <= bottom) && (!done))
    {
      int maxChild = root * 2 + 1;

      if (root*2 == bottom)
	maxChild = root * 2;
      else if (numbers[root * 2] > numbers[root * 2 + 1])
	maxChild = root * 2;
    
      if (numbers[root] < numbers[maxChild])
	{
	  keyType temp = numbers[root];
	  numbers[root] = numbers[maxChild];
	  numbers[maxChild] = temp;
	
	  uint temp2 = data[root];
	  data[root] = data[maxChild];
	  data[maxChild] = temp2;
	
	  root = maxChild;
	}
      else
	done = 1;
    }
}

__kernel void heapSortData(__global keyType* numbers,  __global uint* data, uint array_size)
{
  int i;
 
  for (i = (array_size / 2) - 1; i >= 0; i--)
    siftDownData(numbers, data, i, array_size - 1);
 
  for (i = array_size-1; i >= 1; i--)
    {
      keyType temp = numbers[0];
      numbers[0] = numbers[i];
      numbers[i] = temp;

      uint temp2 = data[0];
      data[0] = data[i];
      data[i] = temp2;

      siftDownData(numbers, data, 0, i-1);
    }
});
      }
    };
  }
}
#undef STRINGIFY
