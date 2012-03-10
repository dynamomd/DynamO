#include <iostream>
#include <vector>
#include <cstdlib>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

int main(int argc, char *argv[]) {
	try {
		// get platforms
		std::vector<cl::Platform> platforms;
		cl::Platform::get(&platforms);

		// Iterate over all platforms
		for(std::vector<cl::Platform>::const_iterator pltfmIt =
			platforms.begin(); pltfmIt != platforms.end(); ++pltfmIt) {
			// print out platform information			
			std::cout << "OpenCL Platform: " << pltfmIt->getInfo<CL_PLATFORM_NAME>() << std::endl
				<< "    Vendor:     " << pltfmIt->getInfo<CL_PLATFORM_VENDOR>() << std::endl
				<< "    Version:    " << pltfmIt->getInfo<CL_PLATFORM_VERSION>() << std::endl
				<< "    Extensions: " << pltfmIt->getInfo<CL_PLATFORM_EXTENSIONS>() << std::endl << std::endl;

			// get platform devices
			std::vector<cl::Device> devices;
			pltfmIt->getDevices(CL_DEVICE_TYPE_ALL, &devices);
			
			// iterate over all devices
			for(std::vector<cl::Device>::const_iterator dvIt = 
				devices.begin(); dvIt != devices.end(); ++dvIt) {
				// print out device information
				std::cout << "    device:   " << dvIt->getInfo<CL_DEVICE_NAME>() << std::endl
				          << "      type:           " << (dvIt->getInfo<CL_DEVICE_TYPE>()==CL_DEVICE_TYPE_CPU?"CPU":"No CPU") << std::endl
				          << "      vendor:         " << dvIt->getInfo<CL_DEVICE_VENDOR>() << std::endl
				          << "      compute units:  " << dvIt->getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << std::endl
				          << "      global memory:  " << dvIt->getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/1048576 << " MB" << std::endl
				          << "      local memory:   " << dvIt->getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/1024 << " KB" << std::endl
				          << "      clock frequency:" << dvIt->getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << " MHz" << std::endl
				          << "      image support:  " << (dvIt->getInfo<CL_DEVICE_IMAGE_SUPPORT>()==CL_TRUE?"Yes":"No") << std::endl
				          << "      extensions:     " << dvIt->getInfo<CL_DEVICE_EXTENSIONS>() << std::endl;
			}
		}
	} catch(cl::Error& err) {
		std::cerr << "OpenCL error: " << err.what() << "(" << err.err() <<
			")" << std::endl;

		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
