#include "lines.hpp"

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include <coil/RenderObj/Arrows.hpp>
# include "../liouvillean/OrientationL.hpp"

magnet::thread::RefPtr<RenderObj>& 
SpLines::getCoilRenderObj() const
{
  M_throw() << "Fix the race condition in here (migrate the particleData and particleColorData to a shared object)";
  if (!_renderObj.isValid())
    {
      _renderObj = new RArrows(range->size(), "Species: " + spName);
      particleData.resize(range->size()*6);
      particleColorData.resize(range->size());
    }

  return _renderObj;
}

void
SpLines::sendRenderData(magnet::CL::CLGLState& CLState) const
{
  CLState.getCommandQueue().enqueueWriteBuffer
    (static_cast<RArrows&>(*_renderObj).getPointData(),
     false, 0, 3 * range->size() * sizeof(cl_float), &particleData[0]);
  
  CLState.getCommandQueue().enqueueWriteBuffer
    (static_cast<RArrows&>(*_renderObj).getDirectionData(),
     false, 0, 3 * range->size() * sizeof(cl_float), &particleData[3*range->size()]);
}

void
SpLines::updateRenderData(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  double diam = getIntPtr()->maxIntDist();
      
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

  _renderObj->getQueue()->queueTask(magnet::function::Task::makeTask(&SpLines::sendRenderData, this, 
								     CLState));
}

void 
SpLines::updateColorObj(magnet::CL::CLGLState&) const
{

}

#endif

