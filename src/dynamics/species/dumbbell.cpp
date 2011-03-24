#ifdef DYNAMO_visualizer
# include "renderobjs/spheres.hpp"
# include <magnet/thread/mutex.hpp>
# include "../liouvillean/OrientationL.hpp"
# include "../interactions/dumbbells.hpp"
# include <magnet/color/HSV.hpp>
# include "dumbbell.hpp"
# include "../liouvillean/CompressionL.hpp"

magnet::thread::RefPtr<RenderObj>& 
SpDumbbells::getCoilRenderObj() const
{
  if (!_renderObj.isValid())
    {
      if (dynamic_cast<const IDumbbells*>(getIntPtr()) == NULL)
	M_throw() << "You must use the IDumbbells interaction for the Dumbbells species type";
      
      _renderObj = new SphereParticleRenderer(2 * range->size(), "Species: " + spName,
					      magnet::function::MakeDelegate(this, &SpDumbbells::updateColorObj),
					      2);
      _coil = new CoilRegister;
    }
  
  return _renderObj;
}

void
SpDumbbells::updateRenderData(magnet::CL::CLGLState& CLState) const
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  std::vector<cl_float4>& particleData 
    = _renderObj.as<SphereParticleRenderer>()._particleData;

  ///////////////////////POSITION DATA UPDATE
  //Check if the system is compressing and adjust the radius scaling factor
  float factor = 1;
  if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
    factor = (1 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime);
 
  double diam = static_cast<const IDumbbells&>(*getIntPtr()).getDiameter();
  double spacing = static_cast<const IDumbbells&>(*getIntPtr()).getLength();
  
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

  if (_renderObj.as<SphereParticleRenderer>().getRecolorOnUpdate())
    updateColorObj(CLState);
  
  _coil->getInstance().getTaskQueue().queueTask(magnet::function::Task::makeTask(&SphereParticleRenderer::sendRenderData, 
										 &(_renderObj.as<SphereParticleRenderer>()), CLState));
}
#endif

