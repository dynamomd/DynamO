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
#include <magnet/CL/scan.hpp>
#include <sstream>
#include <fstream>

#define STRINGIFY(A) #A

namespace magnet {
  namespace CL {
    template<class T>
    class radixSortAMD : public detail::Program
    {
      cl::Kernel _histogramKernel, 
	_permuteKernel,
	_dataPermuteKernel;

      scan<cl_uint> _scanFunctor;
    
      cl::Buffer _buckets, _offsets, _doubleBuffer, _dataDoubleBuffer;
    
      cl_uint _lastSize, _lastDataSize;

    public:
      radixSortAMD():
	_lastSize(0),
	_lastDataSize(0)
      {}

      void build(cl::CommandQueue queue, cl::Context context)
      {
	//Build the scan functor
	_scanFunctor.build(queue, context);

	//And build this functor
	Program::build(queue, context, "");

	_histogramKernel 
	  = cl::Kernel(_program, "histogram");
	_permuteKernel 
	  = cl::Kernel(_program, "permute");
	_dataPermuteKernel 
	  = cl::Kernel(_program, "datapermute");
      }

      void operator()(cl::Buffer keyInput, cl::Buffer keyOutput, cl_uint bits_to_sort = 0)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 64;
	cl_uint bitsPerPass = 4;
	cl_uint maxRadixDigit = 1 << bitsPerPass;
	cl_uint keysPerWorkitem = 256;
	cl_uint nWorkGroups = size / (groupSize * keysPerWorkitem);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % (groupSize * keysPerWorkitem))
	  M_throw() << "Radix sort works on whole multiples of " 
		    << groupSize * keysPerWorkitem << " elements only, please pad your data";

	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * groupSize * maxRadixDigit);
	    
	    _doubleBuffer = cl::Buffer(_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }

	cl::KernelFunctor _histogramFunc
	  = _histogramKernel.bind(_queue, 
				  cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
      
	cl::KernelFunctor _permuteFunc
	  = _permuteKernel.bind(_queue, 
				cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
       
	try {
	  _queue.enqueueCopyBuffer(keyInput, keyOutput, 0, 0, 
									 sizeof(T) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If keyInput and keyOutput are equal no need to copy
	  }

	static int number = 0;
	++number;
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _histogramFunc(keyOutput, 
			   _buckets, 
			   startBit,
			   cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			   size, keysPerWorkitem, bitsPerPass
			   );

	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	    

	    _permuteFunc(keyOutput, _buckets, startBit, 
			 cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			 _doubleBuffer,
			 size, keysPerWorkitem, bitsPerPass
			 );
	    
	    _queue.enqueueCopyBuffer(_doubleBuffer, keyOutput, 0, 0, 
									   sizeof(T) * size);
	  }

	_lastSize = size;
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput,
		      cl::Buffer keyOutput, cl::Buffer dataOutput, 
		      cl_uint bits_to_sort = 0)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 64;
	cl_uint bitsPerPass = 4;
	cl_uint maxRadixDigit = 1 << bitsPerPass;
	cl_uint keysPerWorkitem = 256;
	cl_uint nWorkGroups = size / (groupSize * keysPerWorkitem);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % (groupSize * keysPerWorkitem))
	  M_throw() << "Radix sort works on whole multiples of " 
		    << groupSize * keysPerWorkitem << " elements only, please pad your data";

	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * groupSize * maxRadixDigit);
	    
	    _doubleBuffer = cl::Buffer(_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }

	if (_lastDataSize != size)
	  {
	    _dataDoubleBuffer = cl::Buffer(_context, 
					   CL_MEM_READ_WRITE, sizeof(cl_uint) * size);
	  }

	cl::KernelFunctor _histogramFunc
	  = _histogramKernel.bind(_queue, 
				  cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
      
	cl::KernelFunctor _dataPermuteFunc
	  = _dataPermuteKernel.bind(_queue, 
				    cl::NDRange(nWorkGroups * groupSize), cl::NDRange(groupSize));
       
	try {
	  _queue.enqueueCopyBuffer(keyInput, keyOutput, 0, 0, 
									 sizeof(T) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If keyInput and keyOutput are equal no need to copy
	  }

	try {
	  _queue.enqueueCopyBuffer(dataInput, dataOutput, 0, 0, 
									 sizeof(cl_uint) * size);
	} catch (cl::Error& err)
	  {
	    if (err.err() != CL_MEM_COPY_OVERLAP) throw; //If dataInput and dataOutput are equal no need to copy
	  }

	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _histogramFunc(keyOutput, 
			   _buckets, 
			   startBit,
			   cl::__local(sizeof(cl_ushort) * maxRadixDigit * groupSize),
			   size, keysPerWorkitem, bitsPerPass
			   );

	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	    

	    _dataPermuteFunc(keyOutput, _buckets, dataOutput, startBit, 
			     cl::__local(sizeof(cl_uint) * maxRadixDigit * groupSize),
			     _doubleBuffer, _dataDoubleBuffer,
			     size, keysPerWorkitem, bitsPerPass
			     );
	    
	    _queue.enqueueCopyBuffer(_doubleBuffer, keyOutput, 0, 0, 
									   sizeof(T) * size);
	    _queue.enqueueCopyBuffer(_dataDoubleBuffer, dataOutput, 0, 0, 
									   sizeof(cl_uint) * size);
	  }

	_lastDataSize = _lastSize = size;
      }

      virtual std::string initKernelSrc()
      {
	return
	  "\n#define keyType " 
	  + detail::traits<typename detail::traits<T>::bitshiftable_type>::kernel_type()
	  + "\n#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable\n"
	  STRINGIFY(
/**
 * @brief   Calculates block-histogram bin whose bin size is 256
 * @param   unsortedData    array of unsorted elements
 * @param   buckets         histogram buckets    
 * @param   shiftCount      shift count
 * @param   sharedArray     shared array for thread-histogram bins
  */
__kernel
void histogram(__global const uint* unsortedData,
               __global uint* buckets,
               uint shiftCount,
               __local ushort* sharedArray,
	       uint N,
	       const uint itemsPerThread, const uint radix)
{
//  const uint itemsPerThread = 256;
//  const uint radix = 8;
  const uint radices = 1 << radix;
  const uint mask = radices - 1;

  for (uint blockOffset = get_group_id(0) * get_local_size(0);
       blockOffset * itemsPerThread < N;
       blockOffset += get_global_size(0))
    {
      uint globalID = blockOffset + get_local_id(0);
      uint globalSize = N / itemsPerThread;

      /* Initialize shared array to zero */
      for(int i = 0; i < radices; ++i)
	sharedArray[get_local_id(0) * radices + i] = 0;
      
      barrier(CLK_LOCAL_MEM_FENCE);
      
      /* Calculate thread-histograms */
      for(int i = 0; i < itemsPerThread; ++i)
	{
	  uint value = unsortedData[globalID * itemsPerThread + i] >> shiftCount;
	  value &= mask;
	  sharedArray[get_local_id(0) * radices + value]++;
	}
      
      barrier(CLK_LOCAL_MEM_FENCE);
      
      /* Copy calculated histogram bin to global memory */
      for(int i = 0; i < radices; ++i)
	{
	  uint bucketPos = i * globalSize + globalID;
	  buckets[bucketPos] = sharedArray[get_local_id(0) * radices + i];
	}
    }
}

/**
 * @brief   Permutes the element to appropriate places based on
 *          prescaned buckets values
 * @param   unsortedData        array of unsorted elments
 * @param   scanedBuckets       prescaned buckets for permuations
 * @param   shiftCount          shift count
 * @param   sharedBuckets       shared array for scaned buckets
 * @param   sortedData          array for sorted elements
 */
__kernel
void permute(__global const uint* unsortedKeys,
             __global const uint* scanedBuckets,
             uint shiftCount,
             __local uint* sharedBuckets,
             __global uint* sortedKeys,
	     uint N,
	     const uint itemsPerThread, const uint radix)
{
  const uint radices = 1 << radix;
  const uint mask = radices - 1;
  
  for (uint blockOffset = get_group_id(0) * get_local_size(0);
       blockOffset * itemsPerThread < N;
       blockOffset += get_global_size(0))
    {
      uint globalID = blockOffset + get_local_id(0);
      uint globalSize = N / itemsPerThread;
      
      /* Copy prescaned thread histograms to corresponding thread shared block */
      for(int i = 0; i < radices; ++i)
	{
	  uint bucketPos = i * globalSize + globalID;
	  sharedBuckets[get_local_id(0) * radices + i] = scanedBuckets[bucketPos];
	}

      barrier(CLK_LOCAL_MEM_FENCE);
      
      /* Premute elements to appropriate location */
      for(int i = 0; i < itemsPerThread; ++i)
	{
	  uint value = unsortedKeys[globalID * itemsPerThread + i];
	  value = (value >> shiftCount) & mask;
	  uint index = sharedBuckets[get_local_id(0) * radices + value];
	  sortedKeys[index] = unsortedKeys[globalID * itemsPerThread + i];
	  sharedBuckets[get_local_id(0) * radices + value] = index + 1;
	  barrier(CLK_LOCAL_MEM_FENCE);
	}
    }
}

__kernel
void datapermute(__global const uint* unsortedKeys,
		 __global const uint* scanedBuckets,
		 __global const uint* unsortedData,
		 uint shiftCount,
		 __local uint* sharedBuckets,
		 __global uint* sortedKeys,
		 __global uint* sortedData,
		 uint N,
		 const uint itemsPerThread, const uint radix)
{
  const uint radices = 1 << radix;
  const uint mask = radices - 1;

  for (uint blockOffset = get_group_id(0) * get_local_size(0);
       blockOffset * itemsPerThread < N;
       blockOffset += get_global_size(0))
    {
      uint globalID = blockOffset + get_local_id(0);
      uint globalSize = N / itemsPerThread;

      /* Copy prescaned thread histograms to corresponding thread shared block */
      for(int i = 0; i < radices; ++i)
	{
	  uint bucketPos = i * globalSize + globalID;
	  sharedBuckets[get_local_id(0) * radices + i] = scanedBuckets[bucketPos];
	}

      barrier(CLK_LOCAL_MEM_FENCE);
      
      /* Premute elements to appropriate location */
      for(int i = 0; i < itemsPerThread; ++i)
	{
	  uint value = unsortedKeys[globalID * itemsPerThread + i];
	  value = (value >> shiftCount) & mask;
	  uint index = sharedBuckets[get_local_id(0) * radices + value];
	  sortedKeys[index] = unsortedKeys[globalID * itemsPerThread + i];
	  sortedData[index] = unsortedData[globalID * itemsPerThread + i];
	  sharedBuckets[get_local_id(0) * radices + value] = index + 1;
	  barrier(CLK_LOCAL_MEM_FENCE);
	}
    }
});	
      }
    };
  }
}
#undef STRINGIFY
