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

#include <magnet/cloneptr.hpp>
#include "../../base/is_base.hpp"
#include "../ranges/1range.hpp"
#include "../units/units.hpp"
#include <string>
#include <magnet/thread/refPtr.hpp>

#ifdef DYNAMO_visualizer
# include <coil/RenderObj/RenderObj.hpp>
# include <vector>
#endif

class XMLNode;
namespace xml
{
  class XmlStream;
}
class Particle;
class Interaction;


class Species:public DYNAMO::SimBase
{
public:  
  Species(DYNAMO::SimData*, CRange*, double nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  
  Species(const XMLNode&, DYNAMO::SimData*, unsigned int ID);

  virtual ~Species() {}

  bool isSpecies(const Particle &) const;
  
  const double& getMass() const { return mass; }
  
  unsigned long getCount() const;
  
  unsigned int getID() const { return ID; }
  
  virtual void operator<<(const XMLNode&);

  virtual void initialise();

  friend xml::XmlStream& operator<<(xml::XmlStream&, const Species&);
  
  const std::string& getName() const { return spName; }

  const std::string& getIntName() const { return intName; }

  const Interaction* getIntPtr() const;

  void setIntPtr(Interaction*);

  const magnet::ClonePtr<CRange>& getRange() const { return range; }

  virtual Species* Clone() const { return new Species(*this); }

  virtual double getScalarMomentOfInertia() const { return 0; }

  static Species* getClass(const XMLNode&, DYNAMO::SimData*, unsigned int);

#ifdef DYNAMO_visualizer
  virtual magnet::thread::RefPtr<RenderObj>& getCoilRenderObj() const;

  virtual void updateRenderObj(magnet::CL::CLGLState&) const;

  virtual void updateColorObj(magnet::CL::CLGLState&) const;
#endif

  void setConstantColor(unsigned char r, unsigned char g, unsigned char b)
  { 
    _constColor[0] = r; _constColor[1] = g; _constColor[2] = b; 
    _colorMode = CONSTANT;
  }

protected:
  Species(DYNAMO::SimData*, std::string name, 
	   CRange*, double nMass, std::string nName, 
	   unsigned int ID, std::string nIName="Bulk");
  

  virtual void outputXML(xml::XmlStream&) const;
  
  double mass;

  magnet::ClonePtr<CRange> range;

  std::string spName;
  std::string intName;

  Interaction* IntPtr;

  unsigned int ID;

#ifdef DYNAMO_visualizer
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
};

class SpInertia: public Species
{
public:
  SpInertia(DYNAMO::SimData* sim, std::string name, 
	    CRange* r, double nMass, std::string nName, 
	    unsigned int ID, std::string nIName="Bulk"):
    Species(sim, name, r, nMass, nName, ID, nIName)
  {}

  SpInertia(const XMLNode& XML, DYNAMO::SimData* Sim, unsigned int ID):
    Species(XML,Sim,ID)
  {}
};
