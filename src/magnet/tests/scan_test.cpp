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

#include <magnet/scan.hpp>
#include <iostream>

void runTest(cl::Context context, cl::CommandQueue queue)
{
  magnet::scan<cl_uint> scanFunctor(queue, context);
}

int main(int argc, char *argv[])
{
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
  
  return 0; //Test passed
}

