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
#include <magnet/CL/heapSort.hpp>
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
	std::cout << "Error i = " << i 
		  << " output = " << output[i]
		  << " answer = " << answer[i]
		  << "\n";
	result = false;
      }

  return result;
}

template<class T>
void runTestType(cl::Context context, cl::CommandQueue queue)
{
  cl_uint size = 1024 * 5;

  std::vector<T> input(size);

  std::cout << "##Testing heap sort for " << input.size() << " elements and type " 
	    << magnet::CL::detail::traits<T>::kernel_type()
	    << std::endl;
  
  for(size_t i = 0; i < input.size(); ++i)
    input[i] = input.size() - i - 1;
  
  // create input buffer using pinned memory
  cl::Buffer bufferIn(context, CL_MEM_ALLOC_HOST_PTR |
		      CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE, 
		      sizeof(T) * input.size(), &input[0])
    ;
  
  magnet::CL::heapSort<T> heapSortFunctor(queue, context);

  heapSortFunctor(bufferIn);

  std::vector<T> output(size);
 
  queue.enqueueReadBuffer(bufferIn, CL_TRUE, 0, input.size() *
			  sizeof(T), &output[0]);

  if (!testOutput(input, output))
    M_throw() << "Incorrect output for keysort with size " 
	      << input.size()
	      << " and type "
	      << magnet::CL::detail::traits<T>::kernel_type();

  //Now test with some data!
  //Refresh the input array
  queue.enqueueWriteBuffer(bufferIn, CL_TRUE, 0, input.size() *
			   sizeof(T), &input[0]);

  std::cout << "##key Sort passed, trying data sort\n";

  //Write a data array
  std::vector<cl_uint> data(size);
  for(size_t i = 0; i < input.size(); ++i)
    data[i] = i;

  cl::Buffer dataIn(context, CL_MEM_ALLOC_HOST_PTR |
		    CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE, 
		    sizeof(cl_uint) * data.size(), &data[0])
    ;

  heapSortFunctor(bufferIn, dataIn);
  
  queue.enqueueReadBuffer(dataIn, CL_TRUE, 0, data.size() *
			  sizeof(cl_uint), &data[0]);

  if (!testOutput(input, output))
    M_throw() << "Incorrect output in keys for data and key sort with size " 
	      << input.size()
	      << " and type "
	      << magnet::CL::detail::traits<T>::kernel_type();

  for(size_t i = 0; i < input.size(); ++i)
    if (data[i] != input.size() - 1 - i)
      M_throw() << "Incorrect output in data for data and key sort with size " 
		<< input.size()
		<< " and type "
		<< magnet::CL::detail::traits<T>::kernel_type();  
}

void runTest(cl::Context context, cl::CommandQueue queue)
{
  runTestType<cl_uint>(context, queue);
  runTestType<cl_int>(context, queue);
  runTestType<cl_float>(context, queue);
}

int main(int argc, char *argv[])
{
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
	  if (devIt->getInfo<CL_DEVICE_TYPE>() != CL_DEVICE_TYPE_CPU) continue;

	  std::cout << "#OpenCL CPU device [" << devIt - allDevices.begin() << "]: " <<
	    devIt->getInfo<CL_DEVICE_NAME>() << std::endl;
	  
	  std::vector<cl::Device> devices;
	  devices.push_back(*devIt);
	  
	  cl::Context context(devices);
	  cl::CommandQueue queue(context, devices.front());
	  
	  runTest(context, queue);
	  
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
  
  return 0; //Test passed
}

