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

#define STRINGIFY(A) #A

namespace magnet {
  namespace CL {
    template<class T>
    class radixSortNVIDIA : public detail::Functor
    {
      cl::Kernel _radixSortKernel, 
	_findRadixOffsetsKernel, 
	_reorderKeysKernel,
	_radixSortDataKernel,
	_reorderKeysDataKernel;

      scan<cl_uint> _scanFunctor;
    
      cl::Buffer _buckets, _offsets, _doubleBuffer, _dataDoubleBuffer;
    
      cl_uint _lastSize, _lastDataSize;

      //cl::KernelFunctor _dataSort;

    public:
      radixSortNVIDIA():
	_lastSize(0),
	_lastDataSize(0)
      {}

      void build (cl::CommandQueue queue, cl::Context context)
      {
	//Build the scan functor
	_scanFunctor.build(queue, context);

	//And build this functor
	Functor::build(queue, context, "");

	_radixSortKernel 
	  = cl::Kernel(_program, "radixBlockSortKernel");
	_findRadixOffsetsKernel 
	  = cl::Kernel(_program, "findRadixOffsetsKernel");
	_reorderKeysKernel 
	  = cl::Kernel(_program, "reorderKeys");
	_radixSortDataKernel 
	  = cl::Kernel(_program, "radixBlockSortDataKernel");
	_reorderKeysDataKernel 
	  = cl::Kernel(_program, "reorderKeysData");
      }

      void operator()(cl::Buffer keyInput, cl::Buffer keyOutput, cl_uint bits_to_sort = 0, cl_uint bitsPerPass = 4)
      {
	if (bits_to_sort == 0) 
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);
	cl_uint groupSize = 256;
	cl_uint nWorkGroups = ((size / 4) + groupSize - 1) / groupSize;
	cl_uint maxRadixDigit = 2 << (bitsPerPass - 1);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % 1024)				
	  M_throw() << "Radix sort works on whole multiples of 1024 elements only, please pad your data";

	cl::KernelFunctor _blockSortFunc
	  = _radixSortKernel.bind(_queue, 
				  cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _findRadixOffsetsFunc
	  = _findRadixOffsetsKernel.bind(_queue, 
					 cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _ReorderFunc
	  = _reorderKeysKernel.bind(_queue, 
				    cl::NDRange(size/4), cl::NDRange(groupSize));
       
	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _offsets = cl::Buffer(_context, 
				  CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _doubleBuffer = cl::Buffer(_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }
      
      
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _blockSortFunc(keyInput, _doubleBuffer, size, startBit, bitsPerPass);
	  
	    _findRadixOffsetsFunc(_doubleBuffer, _buckets, _offsets, 
				  size, startBit, bitsPerPass,
				  cl::__local(sizeof(cl_uint) * maxRadixDigit));
	  
	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	  
	    _ReorderFunc(_doubleBuffer, keyOutput, _buckets, _offsets, 
			 size, startBit, bitsPerPass,
			 cl::__local(sizeof(cl_uint) * maxRadixDigit),
			 cl::__local(sizeof(cl_uint) * maxRadixDigit)
			 );
	  }

	_lastSize = size;
      }

      void operator()(cl::Buffer keyInput, cl::Buffer dataInput,
		      cl::Buffer keyOutput, cl::Buffer dataOutput, 
		      cl_uint bits_to_sort = 0, cl_uint bitsPerPass = 4)
      {
	if (bits_to_sort == 0)
	  bits_to_sort = sizeof(T) * 8;

	cl_uint size = keyInput.getInfo<CL_MEM_SIZE>() / sizeof(T);

	if ((dataInput.getInfo<CL_MEM_SIZE>() / sizeof(cl_uint)) != size)
	  M_throw() << "Key and data set size mismatch";

	cl_uint groupSize = 256;
	cl_uint nWorkGroups = ((size / 4) + groupSize - 1) / groupSize;
	cl_uint maxRadixDigit = 2 << (bitsPerPass - 1);

	if (bits_to_sort % bitsPerPass)
	  M_throw() << "The number of bits_to_sort must be a whole multiple of bitsPerPass";

	if (size % 1024)				
	  M_throw() << "Radix sort works on whole multiples of 1024 elements only, please pad your data";

	cl::KernelFunctor _blockSortDataFunc
	  = _radixSortDataKernel.bind(_queue, 
				      cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _findRadixOffsetsFunc
	  = _findRadixOffsetsKernel.bind(_queue, 
					 cl::NDRange(size/4), cl::NDRange(groupSize));
      
	cl::KernelFunctor _reorderDataFunc
	  = _reorderKeysDataKernel.bind(_queue, 
					cl::NDRange(size/4), cl::NDRange(groupSize));
      
	//Create the buffer holding the bit block offsets
	if (_lastSize != size)
	  {
	    _buckets =  cl::Buffer(_context, 
				   CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _offsets = cl::Buffer(_context, 
				  CL_MEM_READ_WRITE, sizeof(cl_uint) * nWorkGroups * maxRadixDigit);
	  
	    _doubleBuffer = cl::Buffer(_context, 
				       CL_MEM_READ_WRITE, sizeof(T) * size);
	  }
 

	if (_lastDataSize != size)
	  {
	    _dataDoubleBuffer = cl::Buffer(_context, 
					   CL_MEM_READ_WRITE, sizeof(cl_uint) * size);
	  }
     
	for (cl_uint startBit = 0; startBit < bits_to_sort; startBit += bitsPerPass)
	  {
	    _blockSortDataFunc(keyInput, dataInput, _doubleBuffer, _dataDoubleBuffer, size, 
			       startBit, bitsPerPass)
	      ;
	  
	    _findRadixOffsetsFunc(_doubleBuffer, _buckets, _offsets, 
				  size, startBit, bitsPerPass,
				  cl::__local(sizeof(cl_uint) * maxRadixDigit))
	      ;
	  
	    //Get the global offsets
	    _scanFunctor(_buckets, _buckets);
	  
	    _reorderDataFunc(_doubleBuffer, _dataDoubleBuffer, keyOutput, dataOutput, 
			     _buckets, _offsets, 
			     size, startBit, bitsPerPass,
			     cl::__local(sizeof(cl_uint) * maxRadixDigit),
			     cl::__local(sizeof(cl_uint) * maxRadixDigit)
			     )
	      ;
	  }

	_lastDataSize = _lastSize = size;
      }

      virtual std::string initKernelSrc()
      {
	return  
	  "#pragma OPENCL EXTENSION cl_khr_byte_addressable_store : enable\n"
	  + scan<cl_uint>().initKernelSrc()
	  + "\n#define keyType4 " 
	  + detail::traits<typename detail::traits<typename detail::traits<T>::bitshiftable_type>::vec4_type>::kernel_type()
	  + "\n#define keyType " 
	  + detail::traits<typename detail::traits<T>::bitshiftable_type>::kernel_type()
	  + "\n"
	  STRINGIFY(
//////////////////////////////////////////////////////////////////////
////////////          Radix Sort Kernels                 /////////////
//////////////////////////////////////////////////////////////////////
//This function will perform a single bit radix sort on a block of data in local memory
void radixBlockSort(keyType4* localkey, uint shift, uint startBit, 
		    __local keyType* share, __local uint* totalTrue)
{
  //First, perform a scan over our little block of 4 items
  uint4 test4;
		  
  test4.x = (((*localkey).x) >> shift) & 0x1;
  test4.y = (((*localkey).y) >> shift) & 0x1;
  test4.z = (((*localkey).z) >> shift) & 0x1;
  test4.w = (((*localkey).w) >> shift) & 0x1;
  
  uint4 sum4 = test4;
  sum4.y += sum4.x;
  sum4.z += sum4.y;
  sum4.w += sum4.z;
  
  __local uint* offsets = (__local uint*)share;
		  
  offsets[get_local_id(0)] = sum4.w;
  
  //Now we must perform a block scan of the sum4.w elements to
  //find the offsets of our little blocks
  scanLocalBlock(offsets, 256, totalTrue);
  
  //We need to know the total number of true statements to calculate offsets for them
  //If something is true, it should go in the second half
  
  sum4 += offsets[get_local_id(0)];
  
  uint totalFalse = 4 * 256 - *totalTrue;
  
  //Now calculate the locations to place each data piece
  sum4.x = (test4.x) 
    ? totalFalse + sum4.x - 1
    : 4 * get_local_id(0) + 0 - sum4.x;

  sum4.y = (test4.y) 
    ? totalFalse + sum4.y - 1
    : 4 * get_local_id(0) + 1 - sum4.y;

  sum4.z = (test4.z) 
    ? totalFalse + sum4.z - 1
    : 4 * get_local_id(0) + 2 - sum4.z;

  sum4.w = (test4.w) 
    ? totalFalse + sum4.w - 1
    : 4 * get_local_id(0) + 3 - sum4.w;
  
  
  //Avoid bank conflicts on the read by interleaving the uint4 data
  sum4.x = (sum4.x & 3) * 256 + (sum4.x >> 2);
  sum4.y = (sum4.y & 3) * 256 + (sum4.y >> 2);
  sum4.z = (sum4.z & 3) * 256 + (sum4.z >> 2);
  sum4.w = (sum4.w & 3) * 256 + (sum4.w >> 2);
  
  //Scatter the data out into a shared memory space      
  
  //AMD needs this extra barrier
  barrier(CLK_LOCAL_MEM_FENCE);
  share[sum4.x] = (*localkey).x;
  share[sum4.y] = (*localkey).y;
  share[sum4.z] = (*localkey).z;
  share[sum4.w] = (*localkey).w;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // The above allows us to read without 4-way bank conflicts:
  (*localkey).x = share[get_local_id(0)];
  (*localkey).y = share[get_local_id(0) +     256];
  (*localkey).z = share[get_local_id(0) + 2 * 256];
  (*localkey).w = share[get_local_id(0) + 3 * 256];
  barrier(CLK_LOCAL_MEM_FENCE);
}

//Just a wrapper kernel function around radixBlockSort, allocates local memory and loads and stores to global 
//memory
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void radixBlockSortKernel(__global keyType4* const keyData, 
			  __global keyType4* outkeyData, 
			  uint array_size, uint startBit, uint nBits)
{
  //Load the 4 keys this thread will sort
  //First load the maximum value
  keyType4 key = keyData[get_global_id(0)];

  //This local array is used twice!
  //First, to perform a scan on 256 uints to determine the
  //little block offsets.
  //Second, to swap the keys between work items 
  __local keyType share[4 * 256];  

  //We have to allocate this inside the kernel!
  __local uint totalTrue;
  
  //Loop over the bits we will sort the block over, and store the totals in the localBuckets array
  for (uint shift = startBit; shift < (startBit + nBits); ++shift)
    radixBlockSort(&key, shift, startBit, share, &totalTrue);
  
  //Write out the sorted keys
  outkeyData[get_global_id(0)] = key;
}

void radixBlockSortData(keyType4* localkey, uint4* localdata, 
			uint shift, uint startBit, 
			__local keyType* share,
			__local uint* datashare,
			__local uint* totalTrue)
{
  //First, perform a scan over our little block of 4 items
  uint4 test4;

  test4.x = (((*localkey).x) >> shift) & 0x1;
  test4.y = (((*localkey).y) >> shift) & 0x1;
  test4.z = (((*localkey).z) >> shift) & 0x1;
  test4.w = (((*localkey).w) >> shift) & 0x1;
  
  uint4 sum4 = test4;
  sum4.y += sum4.x;
  sum4.z += sum4.y;
  sum4.w += sum4.z;
  
  __local uint* offsets = (__local uint*)share;
		  
  offsets[get_local_id(0)] = sum4.w;
  
  //Now we must perform a block scan of the sum4.w elements to
  //find the offsets of our little blocks
  scanLocalBlock(offsets, 256, totalTrue);
  
  //We need to know the total number of true statements to calculate offsets for them
  //If something is true, it should go in the second half
  
  sum4 += offsets[get_local_id(0)];
  
  uint totalFalse = 4 * 256 - *totalTrue;
  
  //Now calculate the locations to place each data piece
  sum4.x = (test4.x) ? totalFalse + sum4.x - 1: 4 * get_local_id(0) + 0 - sum4.x;
  sum4.y = (test4.y) ? totalFalse + sum4.y - 1: 4 * get_local_id(0) + 1 - sum4.y;
  sum4.z = (test4.z) ? totalFalse + sum4.z - 1: 4 * get_local_id(0) + 2 - sum4.z;
  sum4.w = (test4.w) ? totalFalse + sum4.w - 1: 4 * get_local_id(0) + 3 - sum4.w;
  
  
  //Avoid bank conflicts on the read by interleaving the uint4 data
  sum4.x = (sum4.x & 3) * 256 + (sum4.x >> 2);
  sum4.y = (sum4.y & 3) * 256 + (sum4.y >> 2);
  sum4.z = (sum4.z & 3) * 256 + (sum4.z >> 2);
  sum4.w = (sum4.w & 3) * 256 + (sum4.w >> 2);
  
  //Scatter the data out into a shared memory space      
  
  //AMD needs this extra barrier
  barrier(CLK_LOCAL_MEM_FENCE);
  share[sum4.x] = (*localkey).x;
  share[sum4.y] = (*localkey).y;
  share[sum4.z] = (*localkey).z;
  share[sum4.w] = (*localkey).w;
  datashare[sum4.x] = (*localdata).x;
  datashare[sum4.y] = (*localdata).y;
  datashare[sum4.z] = (*localdata).z;
  datashare[sum4.w] = (*localdata).w;
  barrier(CLK_LOCAL_MEM_FENCE);
  
  // The above allows us to read without 4-way bank conflicts:
  (*localkey).x = share[get_local_id(0) + 0 * 256];
  (*localkey).y = share[get_local_id(0) + 1 * 256];
  (*localkey).z = share[get_local_id(0) + 2 * 256];
  (*localkey).w = share[get_local_id(0) + 3 * 256];
  (*localdata).x = datashare[get_local_id(0) + 0 * 256];
  (*localdata).y = datashare[get_local_id(0) + 1 * 256];
  (*localdata).z = datashare[get_local_id(0) + 2 * 256];
  (*localdata).w = datashare[get_local_id(0) + 3 * 256];
  barrier(CLK_LOCAL_MEM_FENCE);
}

//Just a wrapper kernel function around radixBlockSort, allocates local memory and loads and stores to global 
//memory
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void radixBlockSortDataKernel(__global keyType4* const keyData,
			      __global uint4* const dataData,
			      __global keyType4* outkeyData, 
			      __global uint4* outdataData, 
			      uint array_size, uint startBit, uint nBits)
{
  //Load the 4 keys this thread will sort
  //First load the maximum value
  keyType4 key = keyData[get_global_id(0)];
  uint4 data = dataData[get_global_id(0)];

  //This local array is used twice!
  //First, to perform a scan on 256 uints to determine the
  //little block offsets.
  //Second, to swap the keys between work items 
  __local keyType share[4 * 256];  

  __local uint datashare[4 * 256];  

  //We have to allocate this inside the kernel!
  __local uint totalTrue;
  
  //Loop over the bits we will sort the block over, and store the totals in the localBuckets array
  for (uint shift = startBit; shift < (startBit + nBits); ++shift)
    radixBlockSortData(&key, &data, shift, startBit, share, datashare, &totalTrue);
  
  //Write out the sorted keys
  outkeyData[get_global_id(0)] = key;
  outdataData[get_global_id(0)] = data;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
//This kernel takes sorted blocks, offsets, and sizes generated by radixBlockSortKernel.
//It reads and scatters the blocks of data to the appropriate locations in ram
__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void reorderKeys(__global const keyType4* keyData, __global keyType* outkey, 
		 __global uint* buckets, __global uint* offsets,
		 uint array_size, uint startBit, uint nBits,
		 __local uint* localOffsets, __local uint* globalOffsets
		 )
{
  __local keyType4 keyBlocks[256];
  __local keyType* keys = (__local keyType*)keyBlocks;

  //Load the key data
  keyBlocks[get_local_id(0)] = keyData[get_global_id(0)];		  
  barrier(CLK_LOCAL_MEM_FENCE);

  //Load the offset and block offset data
  uint maxRadixDigit = (2 << (nBits - 1));
  uint mask = maxRadixDigit -1;

  if(get_local_id(0) < maxRadixDigit)
    {
      globalOffsets[get_local_id(0)] = buckets[get_local_id(0) * get_num_groups(0) + get_group_id(0)];
      localOffsets[get_local_id(0)] = offsets[get_group_id(0) * maxRadixDigit + get_local_id(0)];
    }

  barrier(CLK_LOCAL_MEM_FENCE);
 
  //Now perform a scatter for all four elements
  for (uint i = 0; i < 4; ++i)
    {
      uint radix = (keys[get_local_id(0) + i * 256] >> startBit) & mask;
      uint globalOffset = globalOffsets[radix] + get_local_id(0) + i * 256 - localOffsets[radix];
      if (globalOffset < array_size)
	outkey[globalOffset] = keys[get_local_id(0) + i * 256];
    }
}

__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void reorderKeysData(__global const keyType4* keyData, __global const uint4* dataData,
		     __global keyType* outkey, __global uint* outdata, 
		     __global uint* buckets, __global uint* offsets,
		     uint array_size, uint startBit, uint nBits,
		     __local uint* localOffsets, __local uint* globalOffsets
		     )
{
  __local keyType4 keyBlocks[256];
  __local keyType* keys = (__local keyType*)keyBlocks;

  __local uint4 dataBlocks[256];
  __local uint* data = (__local uint*)dataBlocks;

  //Load the key and data
  keyBlocks[get_local_id(0)] = keyData[get_global_id(0)];
  dataBlocks[get_local_id(0)] = dataData[get_global_id(0)];		  
  barrier(CLK_LOCAL_MEM_FENCE);

  //Load the offset and block offset data
  uint maxRadixDigit = (2 << (nBits - 1));
  uint mask = maxRadixDigit - 1;
		  
  if(get_local_id(0) < maxRadixDigit)
    {
      globalOffsets[get_local_id(0)] = buckets[get_local_id(0) * get_num_groups(0) + get_group_id(0)];
      localOffsets[get_local_id(0)] = offsets[get_group_id(0) * maxRadixDigit + get_local_id(0)];
    }

  barrier(CLK_LOCAL_MEM_FENCE);
 
  //Now perform a scatter for all four elements
  for (uint i = 0; i < 4; ++i)
    {
      uint radix = (keys[get_local_id(0) + i * 256] >> startBit) & mask;
      uint globalOffset = globalOffsets[radix] + get_local_id(0) + i * 256 - localOffsets[radix];  
      outkey[globalOffset] = keys[get_local_id(0) + i * 256];
      outdata[globalOffset] = data[get_local_id(0) + i * 256];
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
void findRadixOffsets(__global uint* buckets, __global uint* offsets,
		      __local const keyType* block, __local uint* radixOffsets, 
		      uint maxRadixDigit)
{
  //Zero the radix offsets
  if (get_local_id(0) < maxRadixDigit)
    radixOffsets[get_local_id(0)] = 0;

  barrier(CLK_LOCAL_MEM_FENCE);
  
  //We check the first element against the previous key in the shared array
  if ((get_local_id(0) > 0) && (block[get_local_id(0)] != block[get_local_id(0) -1]))
    radixOffsets[block[get_local_id(0)]] = get_local_id(0);
  
  for (size_t i = 1; i < 4; ++i)
    if (block[get_local_id(0) + i * 256] != block[get_local_id(0) + i * 256 - 1])
      radixOffsets[block[get_local_id(0) + i * 256]] = get_local_id(0) + i * 256;
  
  barrier(CLK_LOCAL_MEM_FENCE);

  //Write out the offsets into this blocks radix chunks
  if (get_local_id(0) < maxRadixDigit)
    offsets[get_local_id(0)] = radixOffsets[get_local_id(0)];
  
  barrier(CLK_LOCAL_MEM_FENCE);
  
  //Now calculate the size of each block
  if (get_local_id(0) > 0)
    if (block[get_local_id(0)] != block[get_local_id(0) -1])
      radixOffsets[block[get_local_id(0)-1]] 
	= get_local_id(0) - radixOffsets[block[get_local_id(0)-1]];
  
  for (size_t i = 1; i < 4; ++i)
    if (block[get_local_id(0) + i * 256] != block[get_local_id(0) + i * 256 -1]) 
      radixOffsets[block[get_local_id(0) + i * 256 - 1]] 
	= get_local_id(0) + i * 256 
	- radixOffsets[block[get_local_id(0) + i * 256 - 1]];

  //Now get the length of the last digit
  if (get_local_id(0) == 256 - 1)
    radixOffsets[block[4 * 256 - 1]]
      = 4 * 256 - radixOffsets[block[4 * 256 - 1]];

  barrier(CLK_LOCAL_MEM_FENCE);
    
  //Now write out the block bucket sizes
  if (get_local_id(0) < maxRadixDigit)
    buckets[get_local_id(0) * get_num_groups(0) + get_group_id(0)] 
      = radixOffsets[get_local_id(0)];
}

__kernel __attribute__((reqd_work_group_size(256, 1, 1)))
void findRadixOffsetsKernel(__global keyType4* keyData, __global uint* buckets, 
			    __global uint* offsets,
			    uint array_size, uint startBit, uint nBits, 
			    __local uint* radixOffsets)
{
  __local keyType4 keyradix[256];

  //Offset into the global key data
  keyData += get_group_id(0) * 256;

  //Calculate the max radix digit
  uint maxRadixDigit = (2 << (nBits - 1));
  //Calculate the radix mask
  uint mask = maxRadixDigit - 1;

  //Load the key radix into a local store		  
  keyType4 localKey = keyData[get_local_id(0)];

  keyradix[get_local_id(0)].x = (localKey.x >> startBit) & mask;
  keyradix[get_local_id(0)].y = (localKey.y >> startBit) & mask;
  keyradix[get_local_id(0)].z = (localKey.z >> startBit) & mask;
  keyradix[get_local_id(0)].w = (localKey.w >> startBit) & mask;

  //Wait for the keys load to finish
  barrier(CLK_LOCAL_MEM_FENCE);


  //Offset into the offsets arrays
  offsets += maxRadixDigit * get_group_id(0);
    
  //With all the sorted keys in local memory we can do the offset and bucket count/scan
  findRadixOffsets(buckets, offsets, (__local const keyType*)keyradix, 
		   radixOffsets, maxRadixDigit);
});
      }
    };
  }
}
#undef STRINGIFY
