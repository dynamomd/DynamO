#include "lines.hpp"
#include "../BC/BC.hpp"
#include "../liouvillean/liouvillean.hpp"

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include "renderobjs/lines.hpp"

magnet::thread::RefPtr<coil::RenderObj>& 
SpLines::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    _renderObj = new LineParticleRenderer(range->size(), "Species: " + spName);

  return _renderObj;
}

void
SpLines::updateRenderData(magnet::GL::Context& context) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  ///////////////////////POSITION DATA UPDATE
  //Divide by the maximum box length, to have a natural scale for the visualizer
  const double lengthRescale = 1 / Sim->primaryCellSize.maxElement();

  double diam = getIntPtr()->maxIntDist() * lengthRescale;

  std::vector<cl_float>& particleData 
    = _renderObj.as<LineParticleRenderer>()._particleData;
      
  size_t lineID(0);

  BOOST_FOREACH(unsigned long ID, *range)
    {
      //Position
      Vector pos = Sim->particleList[ID].getPosition();
      Sim->dynamics.BCs().applyBC(pos);
      for (size_t i(0); i < NDIM; ++i)
	particleData[3 * lineID + i] = pos[i] * lengthRescale;
      
      //Vector
      Vector orientation 
	= diam * Sim->dynamics.getLiouvillean()
	.getRotData(Sim->particleList[ID]).orientation;
      for (size_t i(0); i < NDIM; ++i)
	particleData[3 * (range->size() + lineID) + i] = orientation[i];

      ++lineID;
    }

  _coil->getInstance().getTaskQueue().queueTask(magnet::function::Task::makeTask(&LineParticleRenderer::sendRenderData, 
										 &(_renderObj.as<LineParticleRenderer>()), context));
}

#endif

