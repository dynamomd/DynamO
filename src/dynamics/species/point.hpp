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

#include "../../extcode/xmlwriter.hpp"
#include "species.hpp"

class SpPoint: public Species
{
public:  
  SpPoint(DYNAMO::SimData* sim, CRange* r, double nmass, std::string nName,
	  unsigned int ID, std::string nIName="Bulk"):
    Species(sim, "SpPoint", r, nmass, nName, ID, nIName),
    _colorMode(IDHSV),
    _initialColorData(true)
  {}
  
  SpPoint(const XMLNode& XML, DYNAMO::SimData* nSim, unsigned int nID):
    Species(nSim, "", NULL, 0, "", nID,""),
    _colorMode(IDHSV),
    _initialColorData(true)
  { operator<<(XML); }
  
  virtual void initialise();

  virtual void operator<<(const XMLNode& XML);

  virtual Species* Clone() const { return new SpPoint(*this); }

  virtual double getScalarMomentOfInertia() const { return 0; }

#ifdef DYNAMO_visualizer
  virtual magnet::thread::RefPtr<RenderObj>& getCoilRenderObj() const;

  virtual void updateRenderObj(magnet::CL::CLGLState&) const;

  virtual void updateColorObj(magnet::CL::CLGLState&) const;
  
  virtual void showControls(Gtk::ScrolledWindow*);
#endif

  void setConstantColor(unsigned char r, unsigned char g, unsigned char b)
  { 
    _constColor[0] = r; _constColor[1] = g; _constColor[2] = b; 
    _colorMode = CONSTANT;
  }

protected:

#ifdef DYNAMO_visualizer
  void gtkControlCallback();
  mutable magnet::thread::RefPtr<RenderObj> _renderObj;
  mutable std::vector<cl_float4> particleData;
  mutable std::vector<cl_uchar4> particleColorData;
#endif

  typedef enum {
    IDHSV,
    CONSTANT
  } colorMode_t;

  colorMode_t _colorMode;
  unsigned char _constColor[4];
  mutable bool _initialColorData;

  virtual void outputXML(xml::XmlStream& XML) const;
};
