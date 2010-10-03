#include <iostream>

#include <coil/clWindow.hpp>
#include <coil/RenderObj/TestWaves.hpp>
#include <coil/RenderObj/Spheres.hpp>

int main(int argc, char** argv)
{
  try {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    
    cl::Platform clplatform = platforms[0];
    
    CoilMaster::getInstance(argc, argv);
    
    CLGLWindow* CLWindow = new CLGLWindow(1024, 1024,//height, width
					  200, 400,//initPosition (x,y)
					  "GLCLWindow",//title
					  clplatform);
    
    //CLWindow.addRenderObj<RTTestWaves>((size_t)1000, 0.0f);

    std::vector<RTSpheres::SphereDetails> sphereDetailLevels;

    //Shadow testing
//    size_t N = 10 * 1024;
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 2, 1024));
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 1, 1024));
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 0, N-2048));

    //Work computer test render
    size_t N = 1024 * 1000;
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 2, 10));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 1, 1000));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 0, 10000));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::octahedron, 0, 200000));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::tetrahedron, 0, N - 200000 - 10000 - 1000 -10));

    //Home laptop test render
//    size_t N = 10 * 1024;
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 2, 10));
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 1, 100));
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 0, 1000));
//    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::octahedron, 0, N - 1000 - 100 - 10));

    RTSpheres& sphereObject = CLWindow->addRenderObj<RTSpheres>((size_t)N, sphereDetailLevels);

    CoilMaster::getInstance().addWindow(CLWindow);

    //Start the render thread
    CoilMaster::getInstance().bootCoil();

    size_t edit = 0;

    int frameTime = CLWindow->getLastFrameTime();
    while (1)
      {
	//A loop delay function
	//struct timespec duree_nanosleep, duree_out;
	//duree_nanosleep.tv_sec=0;
	//duree_nanosleep.tv_nsec=10000000; //10 ms
	//nanosleep(&duree_nanosleep,&duree_out);
	
	//Now for the update test
	if (frameTime == CLWindow->getLastFrameTime()) continue;
	//The screen was redrawn!
	frameTime = CLWindow->getLastFrameTime();

	//Now try getting access to the sphere position data
	cl_float4* sphereDataPtr = sphereObject.writePositionData(CLWindow->getCommandQueue());	
	++edit;	
	sphereDataPtr[0].w = (edit % 2) ? 0.01 : 0.05;

	//Return it
	sphereObject.returnPositionData(CLWindow->getCommandQueue(), sphereDataPtr);
	
      }

    //Wait for the render thread to exit
    CoilMaster::getInstance().waitForShutdown();
    
  } catch(cl::Error& err) {
    std::cerr << "OpenCL error: " << err.what() << "(" << err.err() 
	      << ")" << std::endl;
    
    return EXIT_FAILURE;
  }
  
  return 0;

  
}

