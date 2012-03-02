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

#include <magnet/CL/detail/extension_wrangler.hpp>
#include <magnet/string/formatcode.hpp>
#include <magnet/string/line_number.hpp>

namespace magnet {
  namespace CL {
    namespace detail {
    
      /* \brief An OpenCL functor object.
       *
       * This class makes it easy to generate OpenCL functors.
       */
      class Program
      {
      public:
	/*! \brief Build the kernel source and store the queue and
	 * context.
	 *
	 * \param buildFlags Any compiler options to pass to the OpenCL compiler.
	 */
	inline void build(::cl::CommandQueue queue, cl::Context context, std::string buildFlags = "")
	{
	  _queue = queue;
	  _context = context;

	  cl::Program::Sources sources;

	  std::string extensions;

	  //Load the double extension if it is present
	  if (magnet::CL::detail::detectExtension(context,"cl_khr_fp64"))
	    extensions += "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n";
	  else if (magnet::CL::detail::detectExtension(context,"cl_amd_fp64"))
	    extensions += "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n";

	  sources.push_back(std::make_pair(extensions.c_str(),
					   extensions.size()));

	  //Now load the kernel source
	  if (_kernelSrc.empty())
	    _kernelSrc = magnet::string::format_code(initKernelSrc());
	  
	  sources.push_back(std::make_pair(_kernelSrc.c_str(),
					   _kernelSrc.size()));
	
	  //Attempt to build the source
	  _program = cl::Program(_context, sources);

	  std::vector<cl::Device> devices = _context.getInfo<CL_CONTEXT_DEVICES>();
	  try {
	    _program.build(devices, buildFlags.c_str());
	  } catch(::cl::Error& err) {
	    for(size_t dev = 0; dev < devices.size(); ++dev) {
	      std::string msg = _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[dev]);
	    
	      if (!msg.size())
		continue;
	    
	      M_throw() << "Compilation failed for device " 
			<< devices[dev].getInfo<CL_DEVICE_NAME>()
			<< "\nBuild Log:\n" << msg
			<< "\nProgram Src:\n" << magnet::string::add_line_numbers(_kernelSrc)
		;
	    }
	  }
	}

	/*! \brief Fetch a kernel object out of the program object.
	 */
	::cl::Kernel operator[](std::string kernelName)
	{ return ::cl::Kernel(_program, kernelName.c_str()); }

	/*! \brief Specifies the initial source of the OpenCL
	 * kernel.
	 *
	 * Every derived \ref Functor class needs to override this and
	 * specify the kernel sources.
	 */
	virtual std::string initKernelSrc() = 0;

      protected:
	::cl::Program _program;
	::cl::CommandQueue _queue;
	::cl::Context _context;
	std::string _kernelSrc;
      };
    }
  }
}
