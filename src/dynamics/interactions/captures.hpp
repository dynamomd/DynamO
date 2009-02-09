/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CICapture_H
#define CICapture_H

#include "interaction.hpp"
#include <set>
#include <vector>
#include <boost/tr1/unordered_set.hpp>

//Expose the boost TR1
namespace std {
  using namespace tr1;
}

class CICapture: public CInteraction
{
public:
  CICapture(DYNAMO::SimData*, C2Range*);

  inline size_t getTotalCaptureCount() const
  { return captureMap.size(); }
  
  bool isCaptured(const CParticle&, const CParticle&) const;

  virtual Iflt getInternalEnergy() const = 0;
  
protected:

  mutable std::unordered_set<std::pair<size_t, size_t> > captureMap;

  bool noXmlLoad;

  virtual bool captureTest(const CParticle&, const CParticle&) const = 0;

  void initCaptureMap();

  void loadCaptureMap(const XMLNode&);

  void outputCaptureMap(xmlw::XmlStream&) const;

  void addToCaptureMap(const CParticle&, const CParticle&) const;
  
  void removeFromCaptureMap(const CParticle&, const CParticle&) const;
  
};

#endif
