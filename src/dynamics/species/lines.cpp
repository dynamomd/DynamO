#include "lines.hpp"

#ifdef dynamo_visualizer
# include <magnet/thread/mutex.hpp>
# include "renderobjs/lines.hpp"
# include "../liouvillean/OrientationL.hpp"

magnet::thread::RefPtr<RenderObj>& 
SpLines::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    _renderObj = new LineParticleRenderer(range->size(), "Species: " + spName);

  return _renderObj;
}

void
SpLines::updateRenderData(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  double diam = getIntPtr()->maxIntDist();

  std::vector<cl_float>& particleData 
    = _renderObj.as<LineParticleRenderer>()._particleData;
      
  size_t lineID(0);

  BOOST_FOREACH(unsigned long ID, *range)
    {
      //Position
      Vector pos = Sim->particleList[ID].getPosition();
      Sim->dynamics.BCs().applyBC(pos);
      for (size_t i(0); i < NDIM; ++i)
	particleData[3 * lineID + i] = pos[i];
      
      //Vector
      Vector orientation 
	= diam * static_cast<const LNOrientation&>(Sim->dynamics.getLiouvillean())
	.getRotData(Sim->particleList[ID]).orientation;
      for (size_t i(0); i < NDIM; ++i)
	particleData[3 * (range->size() + lineID) + i] = orientation[i];

      ++lineID;
    }

  _coil->getInstance().getTaskQueue().queueTask(magnet::function::Task::makeTask(&LineParticleRenderer::sendRenderData, 
										 &(_renderObj.as<LineParticleRenderer>()), CLState));
}

#endif

