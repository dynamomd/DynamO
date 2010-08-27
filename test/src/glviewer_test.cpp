#include "glutMaster.hpp"
#include "glutWindow.hpp"
#include "demoWindow.hpp"
#include <iostream>

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
    
    //CLWindow.addRenderObj<RTTestWaves>((size_t)1000, 0.0f);
  
    CLWindow.addRenderObj<RTSpheres>((size_t)125000,
				     Sphere::icosahedron, (size_t)0,
				     Sphere::icosahedron, (size_t)0,
				     (size_t)1000);
    
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

