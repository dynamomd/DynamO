#include "lines.hpp"

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include <coil/RenderObj/Arrows.hpp>
# include "../liouvillean/OrientationL.hpp"

magnet::thread::RefPtr<RenderObj>& 
SpLines::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      _renderObj = new RLines(range->size());
      particleData.resize(range->size()*6);
    }

  return _renderObj;
}

void
SpLines::updateRenderObj(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  double diam = getIntPtr()->maxIntDist();
      
  size_t lineID(0);
  BOOST_FOREACH(unsigned long ID, *range)
    {
      Vector pos = Sim->particleList[ID].getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      Vector orientation 
	= static_cast<const LNOrientation&>(Sim->dynamics.getLiouvillean())
	.getRotData(Sim->particleList[ID]).orientation;

      Vector a = pos - 0.5 * diam * orientation;
      Vector b = pos + 0.5 * diam * orientation;

      for (size_t i(0); i < NDIM; ++i)
	particleData[6 * lineID + i] = a[i];

      for (size_t i(0); i < NDIM; ++i)
	particleData[6 * lineID + 3 + i] = b[i];

      ++lineID;
    }

  {
    cl::GLBuffer& linedata(static_cast<RArrows&>(*_renderObj).getVertexBuffer());
    linedata.acquire(CLState.getCommandQueue());
    
    CLState.getCommandQueue().enqueueWriteBuffer
      ((cl::Buffer)linedata,
       false, 0, 6 * range->size() * sizeof(cl_float), &particleData[0]);
    
    linedata.release(CLState.getCommandQueue());
  }

}
#endif

