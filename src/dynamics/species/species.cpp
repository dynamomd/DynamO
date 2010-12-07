/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <cstring>
#include "include.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"


Species::Species(DYNAMO::SimData* tmp, CRange* nr, double nMass, 
		   std::string nName, unsigned int nID, std::string nIName):
  SimBase(tmp,"Species", IC_blue),
  mass(nMass),range(nr),spName(nName),intName(nIName),IntPtr(NULL),
  ID(nID)
{}

Species::Species(const XMLNode& XML, DYNAMO::SimData* tmp, unsigned int nID):
  SimBase(tmp,"Species", IC_blue),
  mass(1.0),range(NULL),IntPtr(NULL),
  ID(nID)
{ operator<<(XML); }

Species::Species(DYNAMO::SimData* tmp, std::string name, 
		   CRange* nr, double nMass, std::string nName, 
		   unsigned int nID, std::string nIName):
  SimBase(tmp,name, IC_blue),
  mass(nMass),range(nr),spName(nName),intName(nIName),IntPtr(NULL),
  ID(nID)
{}

const Interaction* 
Species::getIntPtr() const
{ 
#ifdef DYNAMO_DEBUG
  if (IntPtr == NULL)
    M_throw() << "Fetching an unset interaction pointer for a species";
#endif

  return IntPtr; 
}

void
Species::setIntPtr(Interaction* nPtr)
{ IntPtr = nPtr; }

void
Species::initialise()
{ 
  if (IntPtr == NULL)
    M_throw() << "Species missing a matching interaction";
}

xml::XmlStream& operator<<(xml::XmlStream& XML, const Species& g)
{
  g.outputXML(XML);
  return XML;
}

void 
Species::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
    
    try {
      mass = boost::lexical_cast<double>(XML.getAttribute("Mass"))
	* Sim->dynamics.units().unitMass();
      spName = XML.getAttribute("Name");
      intName = XML.getAttribute("IntName");
    } 
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in CSpecies";
      }

}

bool 
Species::isSpecies(const Particle &p1) const
{ return range->isInRange(p1); }

void 
Species::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Mass") 
      << mass / Sim->dynamics.units().unitMass()
      << xml::attr("Name") << spName
      << xml::attr("IntName") << intName
      << xml::attr("Type") << "Point"
      << range;
}

unsigned long 
Species::getCount() const
{
  return range->size();
}

#ifdef DYNAMO_visualizer
# include <magnet/thread/mutex.hpp>
# include "../liouvillean/CompressionL.hpp"

magnet::thread::RefPtr<RenderObj>& 
Species::getCoilRenderObj()
{
  if (!_renderObj.isValid())
    {
      _renderObj = new RTSpheres(range->size());
      particleData.resize(range->size());
    }

  return _renderObj;
}

void
Species::updateRenderObj(magnet::CL::CLGLState& CLState)
{
  if (!_renderObj.isValid())
    M_throw() << "Updating before the render object has been fetched";
  
  //Check if the system is compressing and adjust the radius scaling factor
  float factor = 1;
  if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
    factor = (1 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime);
 
  double sysMass = 0;
  Vector COM;

  double diam = getIntPtr()->hardCoreDiam() * factor;
  sysMass += range->size() * getMass();
      
  BOOST_FOREACH(unsigned long ID, *range)
    {
      Vector pos = Sim->particleList[ID].getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      for (size_t i(0); i < NDIM; ++i)
	particleData[ID].s[i] = pos[i];
      
      particleData[ID].w = diam * 0.5;
    }

  {
    CLState.getCommandQueue().enqueueWriteBuffer
      (static_cast<RTSpheres&>(*_renderObj).getSphereDataBuffer(),
       false, 0, range->size() * sizeof(cl_float4), &particleData[0]);
  }

}
#endif

Species* 
Species::getClass(const XMLNode& XML, DYNAMO::SimData* tmp, unsigned int nID)
{
  if (!XML.isAttributeSet("Type"))
    return new Species(XML, tmp, nID);

  if (!std::strcmp(XML.getAttribute("Type"), "Point"))
    return new Species(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "SphericalTop"))
    return new SpSphericalTop(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "Lines"))
    return new SpLines(XML, tmp, nID);
  else if (!std::strcmp(XML.getAttribute("Type"), "FixedCollider"))
    return new SpFixedCollider(XML, tmp, nID);
  else 
    M_throw() << XML.getAttribute("Type")
	      << ", Unknown type of species encountered";

}
