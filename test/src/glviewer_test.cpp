#include <iostream>

#include "clWindow.hpp"
#include "RenderObj/TestWaves.hpp"
#include "RenderObj/Spheres.hpp"

int main(int argc, char** argv)
{  
  try {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    
    cl::Platform clplatform = platforms[0];
    
    GlutMaster glutMaster(argc, argv);
    
    CLGLWindow CLWindow(glutMaster,
			500, 500,//height, width
			200, 400,//initPosition (x,y)
			"GLCLWindow",//title
			clplatform);
    
    CLWindow.addRenderObj<RTTestWaves>((size_t)1000, 0.0f);
  
    std::vector<RTSpheres::SphereDetails> sphereDetailLevels;

    size_t N = 125000;
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 1, 1000));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 0, 10000));
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::tetrahedron, 0, 50000));

    CLWindow.addRenderObj<RTSpheres>((size_t)125000, sphereDetailLevels);
    
    bool noFPSLimit = true;

    int oldTime = glutGet(GLUT_ELAPSED_TIME);
    while (1) { 
      glutMainLoopEvent();
      
      int currTime = glutGet(GLUT_ELAPSED_TIME);
      
      if (noFPSLimit || (currTime - oldTime > 32))
	{
	  GlutMaster::CallBackIdleFunc(); 
	  oldTime = currTime;
	}      
    }
    
  } catch(cl::Error& err) {
    std::cerr << "OpenCL error: " << err.what() << "(" << err.err() 
	      << ")" << std::endl;
    
    return EXIT_FAILURE;
  }
  
  return 0;

  
}

