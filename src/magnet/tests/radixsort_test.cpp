/*    DYNAMO:- Event driven molecular dynamics simulator 
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

#include <iostream>
#include <magnet/CL/radixsort.hpp>
#include <algorithm>

template<class T>
bool testOutput(const std::vector<T>& input, const std::vector<T>& output)
{
  std::vector<T> answer(input);
  
  std::sort(answer.begin(), answer.end());
  
  bool result = true;
  for (size_t i(0); i < output.size(); ++i)
    if (output[i] != answer[i])
      {
	//std::cout << "Error i = " << i 
	//	  << " output = " << output[i]
	//	  << " answer = " << answer[i]
	//	  << "\n";
	result = false;
      }

  return result;
}

template<class T>
bool runTestType(cl::Context context, cl::CommandQueue queue)
{
  cl_uint size = 1024 * 5;

  std::vector<T> input(size);

  std::cout << "##Testing radix sort for " << input.size() << " elements and type " 
	    << magnet::CL::detail::traits<T>::kernel_type();
  
  for(size_t i = 0; i < input.size(); ++i)
    input[i] = input.size() - i - 1;
  
  // create input buffer using pinned memory
  cl::Buffer bufferIn(context, CL_MEM_ALLOC_HOST_PTR |
		      CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE, 
		      sizeof(T) * input.size(), &input[0])
    ;
  
  magnet::CL::radixSort<T> radixSortFunctor;
  radixSortFunctor.build(queue, context);
  radixSortFunctor(bufferIn, bufferIn);

  std::vector<T> output(size);
 
  queue.enqueueReadBuffer(bufferIn, CL_TRUE, 0, input.size() *
			  sizeof(T), &output[0]);

  bool failed = !testOutput(input, output);

  std::cout << " key(only) " << (failed ? "FAILED" : "PASSED") << ", "; 

  //Now test with some data!
  //Refresh the input array
  queue.enqueueWriteBuffer(bufferIn, CL_TRUE, 0, input.size() *
			   sizeof(T), &input[0]);

  //Write a data array
  std::vector<cl_uint> data(size);
  for(size_t i = 0; i < input.size(); ++i)
    data[i] = i;

  cl::Buffer dataIn(context, CL_MEM_ALLOC_HOST_PTR |
		    CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE, 
		    sizeof(cl_uint) * data.size(), &data[0])
    ;

  radixSortFunctor(bufferIn, dataIn, bufferIn, dataIn);
  
  queue.enqueueReadBuffer(dataIn, CL_TRUE, 0, data.size() *
			  sizeof(cl_uint), &data[0]);

  bool keyfail = !testOutput(input, output);

  std::cout << " key " << (keyfail ? "FAILED" : "PASSED"); 

  bool datafail = false;
  for(size_t i = 0; i < input.size(); ++i)
    if (data[i] != input.size() - 1 - i)
      datafail = true;

  std::cout << " data " << (datafail ? "FAILED" : "PASSED") << std::endl;
  
  return failed || keyfail || datafail;
}

bool runTest(cl::Context context, cl::CommandQueue queue)
{
  bool fail = runTestType<cl_uint>(context, queue)
    || runTestType<cl_int>(context, queue)
    || runTestType<cl_float>(context, queue);

  //Special case for the AMD CPU implementation
  //if (magnet::CL::detail::detectExtension(context, "cl_amd_fp64")
  //    && (queue.getInfo<CL_QUEUE_DEVICE>().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU))
  //  {
  //    if (runTestType<cl_double>(context, queue))
  //	std::cerr << "EXPECTED FAIL: AMD double precision radix sort fails due to a bug in the stream SDK.\n";
  //    else
  //	M_throw() << "The bug in the ATI stream SDK has been fixed!";
  //  }
  //else
  //fail = fail || runTestType<cl_double>(context, queue);

  return fail;
}

int main(int argc, char *argv[])
{
  bool fail = false;
  try {
    //Test all devices and platforms for compatability
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    for(std::vector<cl::Platform>::const_iterator pltfmIt =
	  platforms.begin(); pltfmIt != platforms.end(); ++pltfmIt) 
      {
	std::cout << "OpenCL platform [" << pltfmIt - platforms.begin() << "]: " <<
	  pltfmIt->getInfo<CL_PLATFORM_NAME>() << std::endl;
      
	std::vector<cl::Device> allDevices;
	pltfmIt->getDevices(CL_DEVICE_TYPE_ALL, &allDevices);
      
	for(std::vector<cl::Device>::const_iterator devIt = allDevices.begin(); 
	    devIt != allDevices.end(); ++devIt)
	  {
	    std::cout << "#OpenCL device [" << devIt - allDevices.begin() << "]: " <<
	      devIt->getInfo<CL_DEVICE_NAME>() << std::endl;
	  
	    std::vector<cl::Device> devices;
	    devices.push_back(*devIt);
	  
	    cl::Context context(devices);
	    cl::CommandQueue queue(context, devices.front());
	  
	    fail |= runTest(context, queue);
	  
	  }
      }
  } catch(magnet::exception& err) {
    std::cerr << "Magnet error: " << err.what()
	      << std::endl;

    return 1;
    
  } catch(cl::Error& err) {
    std::cerr << "OpenCL error: " << err.what() << "(" << err.err() 
	      << ")" << std::endl;

    return 1;//Test failed
  } 
  
  return fail;
}

