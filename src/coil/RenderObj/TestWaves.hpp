#pragma once
#include "Triangles.hpp"
#include <time.h>

#include "../clWindow.hpp"

class RTTestWaves : public RTriangles
{
public:
  RTTestWaves(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers, size_t N, float Yoffset);

  virtual void clTick(cl::CommandQueue& CmdQ, cl::Context& Context);

protected:
  static const std::string kernelsrc;

  cl::Kernel kernel;
  timespec startTime;

  size_t _N;
  float _Yoffset;
};

template<>
inline void CLGLWindow::addRenderObj<RTTestWaves,size_t,float>(size_t N, float Yoffset)
{
  RenderObjects.push_back(new RTTestWaves(_clcmdq, _clcontext, _cldevice, _hostTransfers, N, Yoffset));
}

