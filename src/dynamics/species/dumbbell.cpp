#ifdef DYNAMO_visualizer
# include "renderobjs/spheres.hpp"
# include <magnet/thread/mutex.hpp>
# include "../liouvillean/OrientationL.hpp"
# include <magnet/color/HSV.hpp>
# include "dumbbell.hpp"


magnet::thread::RefPtr<RenderObj>& 
SpDumbbells::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      _renderObj = new SphereParticleRenderer(2 * range->size(), "Species: " + spName,
					      magnet::function::MakeDelegate(this, &SpDumbbells::updateColorObj));
      particleData.resize(2 * range->size());
      particleColorData.resize(range->size()); //We just queue two copies
    }

  return _renderObj;
}

void 
SpDumbbells::updateColorObj(magnet::CL::CLGLState& CLState) const
{
  SpPoint::updateColorObj(CLState);
  
  {//A second copy to just duplicate the data for the two spheres of the dumbbells
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getColorDataBuffer(),
       false, range->size() * sizeof(cl_uchar4), range->size() * sizeof(cl_uchar4), &particleColorData[0]);
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
	particleData[sphID].s[i] = pos[i];

      pos = cpos - orientation;
      for (size_t i(0); i < NDIM; ++i)
	particleData[range->size() + sphID].s[i] = pos[i];
      
      particleData[sphID].w = diam * 0.5;
      particleData[range->size() + sphID].w = diam * 0.5;
      ++sphID;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getSphereDataBuffer(),
       false, 0, 2 * range->size() * sizeof(cl_float4), &particleData[0]);
  }
}
#endif

