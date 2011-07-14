/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.marcusbannerman.co.uk/dynamo
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
    class scan : public detail::Functor
    {
      cl::Kernel _prescanKernel, _uniformAddKernel;
      std::vector<cl::Buffer> _partialSumBufferStack;
      cl_uint _lastSize;
      
    public:
      scan(): _lastSize(0) {}

      void build(::cl::CommandQueue queue,cl::Context context)
      {
	Functor::build(queue, context, "");
	
	_prescanKernel = cl::Kernel(_program, "prescan");
	_uniformAddKernel = cl::Kernel(_program, "uniformAdd");
      }

      inline void operator()(cl::Buffer input, cl::Buffer output)
      {
	cl_uint size = input.getInfo<CL_MEM_SIZE>() / sizeof(T);

	//Workgroups of 256 work-items process 512 elements
	cl_uint nGroups = (((size / 2) + 256 - 1) / 256);

	//Check if the size has changed from last call,
	//We cache the partial sum buffers as reallocating buffers is slow
	if (size != _lastSize) 
	  {
	    //Rebuild the stack of buffers
	    _partialSumBufferStack.clear();

	    for (cl_uint stageSize = nGroups; stageSize > 1; 
		 stageSize = (stageSize + 511) / 512)
	      _partialSumBufferStack.push_back(cl::Buffer(_context,
							  CL_MEM_READ_WRITE, sizeof(cl_uint) * stageSize));	  

	    _partialSumBufferStack.push_back(cl::Buffer(_context,
							CL_MEM_READ_WRITE, sizeof(cl_uint)));	  
	  }

	//Start the recursion
	recursionFunction(input, output, size, 0);

	_lastSize = size;
      }

      
      virtual std::string initKernelSrc()
      {
	return std::string("#define scantype ") + detail::traits<T>::kernel_type() + "\n"
	  STRINGIFY(
//////////////////////////////////////////////////////////////////////
////////////          Scan Kernels                       /////////////
//////////////////////////////////////////////////////////////////////
		
//Performs a scan on a continuous "block" of "blocksize" elements in local memory
//This function will synchronize the threads before and after the sum
//Also returns the totalSum in a local variable
void scanLocalBlock(__local scantype* block, uint blocksize, __local uint* totalSum)
{
  int stride = 1;
  
  for (int d = blocksize / 2; d > 0; d >>= 1) // build sum in place up the tree
    {
      barrier(CLK_LOCAL_MEM_FENCE);
      
      if (get_local_id(0) < d)
	{
	  int ai = stride*(2*get_local_id(0)+1)-1;
	  block[ai + stride] += block[ai];
	}
      stride *= 2;
    }
		       
  //Store the last element in the totalSum variable
  //and clear the last element.
  //We sync as we're about to read back from local memory, but a judicious
  //choice of the work_item might make this unnecessary
  barrier(CLK_LOCAL_MEM_FENCE);
		       
  if (get_local_id(0) == 0)
    {
      totalSum[0] = block[blocksize - 1];
      block[blocksize - 1] = 0; 
    }
		       
  for (int d = 1; d < blocksize; d *= 2) // traverse down tree & build scan
    {
      stride >>= 1;
      barrier(CLK_LOCAL_MEM_FENCE);
      if (get_local_id(0) < d)
	{
	  int ai = stride*(2*get_local_id(0)+1)-1;
	  scantype t = block[ai];
	  block[ai] = block[ai + stride];
	  block[ai+stride] += t;
	}
    }
		       
  barrier(CLK_LOCAL_MEM_FENCE);
}
		     
//This kernel performs scan's on blocks of (2 * 256) elements.
//The code in the kernel just loads the local array with data, calls
//scanLocalBlock(__local scantype* block), and returns the total sum
//in partial_sums.
		     
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void prescan(__global scantype *g_idata, __global scantype *g_odata, 
	     __global scantype *partial_sums, uint n)
{
  __local scantype localBlock[2 * 256]; 
		       
  size_t memoffset = 2 * 256 * get_group_id(0);
		       
  //Do the offsets
  g_idata += memoffset;
  g_odata += memoffset;
		       
  localBlock[2 * get_local_id(0) + 0] = 0;
  localBlock[2 * get_local_id(0) + 1] = 0;

  //Seems redundant but around 10.11 the ati drivers
  //started misbehaving without it
  barrier(CLK_LOCAL_MEM_FENCE); 
		       
  // load input into shared memory
  if ((memoffset + 2*get_local_id(0)) < n)
    localBlock[2 * get_local_id(0)] =  g_idata[2 * get_local_id(0)]; 
		       
  if ((memoffset + 2 * get_local_id(0) + 1) < n)
    localBlock[2 * get_local_id(0)+1] = g_idata[2 * get_local_id(0) + 1];
		       
  __local uint totalSum;
  scanLocalBlock(localBlock, 2 * 256, &totalSum);
		       
  if (get_local_id(0)==0)
    partial_sums[get_group_id(0)] = totalSum;
		       
  if ((memoffset + 2 * get_local_id(0)) < n)
    g_odata[2*get_local_id(0)] = localBlock[2*get_local_id(0)]; // write results to device memory
		       
  if ((memoffset + 2 * get_local_id(0) + 1) < n)
    g_odata[2 * get_local_id(0) + 1] = localBlock[2*get_local_id(0)+1];
}
		     
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void uniformAdd(__global scantype *g_idata, __global scantype *g_odata, 
		__global scantype *partial_sums, uint n)
{
  if (get_group_id(0) == 0) return;
		       
  __local scantype increment;
		       
  if (get_local_id(0) == 0)
    increment = partial_sums[get_group_id(0)];

  barrier(CLK_LOCAL_MEM_FENCE);
		       
  uint offset = 2 * 256 * get_group_id(0);
  g_odata += offset;
  g_idata += offset;
		       
  if (offset + get_local_id(0) < n)
    g_odata[get_local_id(0)] = g_idata[get_local_id(0)] + increment;
		       
  if (offset + 256 + get_local_id(0) < n)
    g_odata[256 + get_local_id(0)] 
      = g_idata[256 + get_local_id(0)] + increment; 
});
      }

    protected:

      inline void recursionFunction(cl::Buffer input, cl::Buffer output, cl_uint size, cl_uint stage)
      {
	cl_uint nGroups = (((size / 2) + 256 - 1) / 256);

	_prescanKernel.bind(_queue, 
			    cl::NDRange(256 * nGroups), 
			    cl::NDRange(256))
	  (input, output, _partialSumBufferStack[stage], size);
      
      
	//Recurse if we've got more than one workgroup
	if (nGroups > 1)
	  {
	    recursionFunction(_partialSumBufferStack[stage], _partialSumBufferStack[stage], 
			      nGroups, stage+1);
	  
	    _uniformAddKernel.bind(_queue,
				   cl::NDRange(nGroups * 256), 
				   cl::NDRange(256))
	      (output, output, _partialSumBufferStack[stage], size);
	  }
      }
    };
  }
}

#undef STRINGIFY
