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

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>

namespace magnet {
  namespace CL {
    namespace detail {
      //Check if a device supports an OpenCL extension
      inline bool detectExtension(cl::Device device, const std::string extension)
      {
	std::string extensions = device.getInfo<CL_DEVICE_EXTENSIONS>();
	return extensions.find(extension) != std::string::npos;
      }

      //Check if all devices in a context support an OpenCL extension
      inline bool detectExtension(cl::Context context, const std::string extension)
      {
	std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	for (std::vector<cl::Device>::const_iterator devPtr = devices.begin();
	     devPtr != devices.end(); ++devPtr)
	  if (!detectExtension(*devPtr, extension)) return false;

	return true;
      }
    }
  }
}
