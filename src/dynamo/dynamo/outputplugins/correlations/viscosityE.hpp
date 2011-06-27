/*  dynamo:- Event driven molecular dynamics simulator 
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
#include "../outputplugin.hpp"
#include "../../datatypes/vector.hpp"
#include <boost/circular_buffer.hpp>

class OPViscosityE: public OutputPlugin
{
  typedef boost::array<double, NDIM> col;
  typedef boost::array<col, NDIM> matrix;
  
public:
  OPViscosityE(const dynamo::SimData*, const magnet::xml::Node& XML);

  virtual void initialise();

  virtual void output(magnet::xml::XmlStream &);
  
  virtual OutputPlugin* Clone() const { return new OPViscosityE(*this); }
  
  virtual void eventUpdate(const GlobalEvent&, const NEventData&);

  virtual void eventUpdate(const LocalEvent&, const NEventData&);
  
  virtual void eventUpdate(const System&, const NEventData&, const double&);
  
  virtual void eventUpdate(const IntEvent&, const PairEventData&);

  void stream(const double&);

  virtual void operator<<(const magnet::xml::Node&);

protected:
  void impulseDelG(const PairEventData&);
  void impulseDelG(const NEventData&);

  void updateConstDelG(const PairEventData&);
  void updateConstDelG(const ParticleEventData&);
  void updateConstDelG(const NEventData&);
  
  void newG(const matrix&);
  void accPass();

  matrix avgTrace;

  size_t count;
  double dt, currentdt;
  matrix constDelG, delG;

  size_t currlen;
  bool notReady;

  size_t CorrelatorLength;

  boost::circular_buffer<matrix> G;
  std::vector<matrix> accG2;
  double dtfactor;
};
