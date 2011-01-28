#include "dumbbell.hpp"

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include <coil/RenderObj/Spheres.hpp>
# include "../liouvillean/OrientationL.hpp"
# include <magnet/HSV.hpp>

magnet::thread::RefPtr<RenderObj>& 
SpDumbbells::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      _renderObj = new RTSpheres(2 * range->size(), "Species: " + spName);
      particleData.resize(2 * range->size());
      particleColorData.resize(2 * range->size());
    }

  return _renderObj;
}

void 
SpDumbbells::updateColorObj(magnet::CL::CLGLState& CLState) const
{

  switch (_colorMode)
    {
    case IDHSV:
      {
	size_t np = range->size();
	for (size_t sphID(0); sphID < np; ++ sphID)
	  {
	    magnet::color::HSVtoRGB(particleColorData[2 * sphID + 0], ((float)(sphID)) / np);
	    magnet::color::HSVtoRGB(particleColorData[2 * sphID + 1], ((float)(sphID)) / np);
	  }
      }
      break;
    case CONSTANT:
      for (size_t id(0); id < 2 * range->size(); ++id)
	for (size_t cc(0); cc < 4; ++cc)
	  particleColorData[id].s[cc] = _constColor[cc];
      break;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getColorDataBuffer(),
       false, 0, 2 * range->size() * sizeof(cl_uchar4), &particleColorData[0]);
  }
}


void
SpDumbbells::updateRenderObj(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  double diam = getIntPtr()->maxIntDist();
  double spacing = diam;
  
  size_t sphID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      Vector cpos = Sim->particleList[ID].getPosition();
      Vector orientation 
	= 0.5 * spacing * static_cast<const LNOrientation&>(Sim->dynamics.getLiouvillean())
	.getRotData(Sim->particleList[ID]).orientation;
      
      Sim->dynamics.BCs().applyBC(cpos);
      
      Vector pos = cpos + orientation;
      for (size_t i(0); i < NDIM; ++i)
	particleData[2 * sphID + 0].s[i] = pos[i];

      pos = cpos - orientation;
      for (size_t i(0); i < NDIM; ++i)
	particleData[2 * sphID + 1].s[i] = pos[i];
      
      particleData[2 * sphID + 0].w = diam * 0.5;
      particleData[2 * sphID + 1].w = diam * 0.5;
      ++sphID;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getSphereDataBuffer(),
       false, 0, 2 * range->size() * sizeof(cl_float4), &particleData[0]);
  }
}
#endif

