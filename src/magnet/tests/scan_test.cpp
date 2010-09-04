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
#include <magnet/scan.hpp>
#include <numeric>

template<class T>
bool testOutput(const std::vector<T>& input, const std::vector<T>& output)
{
  std::vector<T> answer;
  answer.resize(input.size());
  
  //This is an inclusive scan!
  std::partial_sum(input.begin(), input.end(), answer.begin());
  
  bool result = true;
  for (size_t i(1); i < output.size(); ++i)
    if (output[i] != answer[i-1])
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
  cl_uint size = 1024 * 2 + 512;

  std::vector<T> input(size);

  std::cout << "Testing scan for " << input.size() << " elements and type " 
	    << magnet::detail::traits<T>::kernel_type()
	    << std::endl;
  
  for(size_t i = 0; i < input.size(); ++i)
    input[i] = i+1;
  
  // create input buffer using pinned memory
  cl::Buffer bufferIn(context, CL_MEM_ALLOC_HOST_PTR |
		      CL_MEM_COPY_HOST_PTR | CL_MEM_READ_WRITE, 
		      sizeof(cl_uint) * input.size(), &input[0])
    ;
  
  magnet::scan<T> scanFunctor(queue, context);

  scanFunctor(bufferIn, bufferIn, size);

  std::vector<T> output(size);
 
  queue.enqueueReadBuffer(bufferIn, CL_TRUE, 0, input.size() *
			  sizeof(T), &output[0]);

  if (!testOutput(input, output))
    M_throw() << "Incorrect output for size " 
	      << input.size()
	      << " and type "
	      << magnet::detail::traits<T>::kernel_type();
}

void runTest(cl::Context context, cl::CommandQueue queue)
{
  runTestType<cl_uint>(context, queue);
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
	  std::cout << "##OpenCL device [" << devIt - allDevices.begin() << "]: " <<
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

