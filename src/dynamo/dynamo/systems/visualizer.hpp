/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
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
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <coil/clWindow.hpp>
#include <dynamo/systems/system.hpp>

namespace coil {
class DataSet;
}
namespace dynamo {
class SVisualizer : public System {
public:
  SVisualizer(dynamo::Simulation *, std::string, double);

  virtual NEventData runEvent();

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node &) {}

  void particlesUpdated(const NEventData &);

protected:
  SVisualizer(
      const SVisualizer &); // Cannot copy due to the coil update connection

  virtual void outputXML(magnet::xml::XmlStream &) const {}

  shared_ptr<coil::CLGLWindow> _window;
  coil::CoilRegister _coil;

  void initDataSet();
  void updateRenderData();

  shared_ptr<coil::DataSet> _particleData;
  boost::posix_time::ptime _lastUpdate;
  std::vector<std::vector<GLuint>> _interactionIDs;
};
} // namespace dynamo
#endif
