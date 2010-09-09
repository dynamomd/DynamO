#  Makefile for glutMaster demo
#  version 0.3  
#  George Stetten and Korin Crawford

FILES   = glutMaster.c++ main.c++ glutWindow.c++ demoWindow.c++

OBJECTS = glutMaster.o main.o glutWindow.o demoWindow.o clWindow.o RenderObj/RenderObj.o RenderObj/Triangles.o RenderObj/TestWaves.o RenderObj/Spheres.o RenderObj/Primatives/Sphere.o

EXECFIL = run
LIBS = -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm -lOpenCL -lrt -lGLEW

CXXFLAGS = -O3 -g -Wall -Werror -ansi -I include/

$(EXECFIL): $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) -o $(EXECFIL) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :	
	rm -f *.o *~ run RenderObj/*.o RenderObj/Primatives/*.o
