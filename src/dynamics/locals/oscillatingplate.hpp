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

#include "local.hpp"
#include "../coilRenderObj.hpp"

class CLOscillatingPlate: public Local, public CoilRenderObj
{
public:
  CLOscillatingPlate(const magnet::xml::Node&, DYNAMO::SimData*);
  CLOscillatingPlate(DYNAMO::SimData*, Vector, Vector, double, 
		     double, double, double, double, std::string, CRange*, 
		     double timeshift = 0, bool nstrongPlate = false);

  virtual ~CLOscillatingPlate() {}

  virtual Local* Clone() const { return new CLOscillatingPlate(*this); };

  virtual LocalEvent getEvent(const Particle&) const;

  virtual void runEvent(const Particle&, const LocalEvent&) const;
  
  virtual bool isInCell(const Vector &, const Vector &) const;

  virtual void initialise(size_t);

  virtual void operator<<(const magnet::xml::Node&);

  virtual void write_povray_info(std::ostream&) const;

  Vector getPosition() const;

  Vector getVelocity() const;

  double getPlateEnergy() const;

  const Vector& getCentre() const { return rw0; }

#ifdef DYNAMO_visualizer
  virtual magnet::thread::RefPtr<RenderObj>& getCoilRenderObj() const;
  virtual void updateRenderData(magnet::CL::CLGLState&) const;
#endif

protected:
#ifdef DYNAMO_visualizer
  mutable magnet::thread::RefPtr<RenderObj> _renderObj;
#endif

  virtual void outputXML(xml::XmlStream&) const;

  bool strongPlate;
  Vector rw0;
  Vector nhat;
  double omega0;
  double sigma;
  double e;
  mutable double delta;
  double mass;
  mutable double timeshift;
  mutable size_t lastID;
  mutable long double lastdSysTime;
};
