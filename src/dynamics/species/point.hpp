/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#pragma once

#ifdef DYNAMO_visualizer
# include <coil/coilMaster.hpp>
#endif

#include "species.hpp"
#include <magnet/xmlwriter.hpp>
#include <memory>

namespace Gtk {
  class VBox;
  class RadioButton;
}

class SpPoint: public Species
{
public:  
  SpPoint(DYNAMO::SimData* sim, CRange* r, double nmass, std::string nName,
	  unsigned int ID, std::string nIName="Bulk"):
    Species(sim, "SpPoint", r, nmass, nName, ID, nIName)
  {}
  
  SpPoint(const magnet::xml::Node& XML, DYNAMO::SimData* nSim, unsigned int nID):
    Species(nSim, "", NULL, 0, "", nID,"")
  { operator<<(XML); }

  virtual void initialise();

  virtual void operator<<(const magnet::xml::Node& XML);

  virtual Species* Clone() const { return new SpPoint(*this); }

  virtual double getScalarMomentOfInertia() const { return 0; }

#ifdef DYNAMO_visualizer
  virtual magnet::thread::RefPtr<RenderObj>& getCoilRenderObj() const;
  virtual void updateRenderData(magnet::CL::CLGLState&) const;
  virtual void updateColorObj(magnet::CL::CLGLState&) const;
#endif

protected:

#ifdef DYNAMO_visualizer
  mutable magnet::thread::RefPtr<RenderObj> _renderObj;
  mutable magnet::thread::RefPtr<CoilRegister> _coil;
#endif

protected:
  virtual void outputXML(xml::XmlStream& XML) const;
};
