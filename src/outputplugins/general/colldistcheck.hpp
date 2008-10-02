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

#ifndef COPCollDistCheck_HPP
#define COPCollDistCheck_HPP

#include "../outputplugin.hpp"
#include <map>
#include "../../datatypes/histogram.hpp"

class COPCollMatrix;

class COPCollDistCheck: public COutputPlugin
{
public:
  COPCollDistCheck(const DYNAMO::SimData*);

  ~COPCollDistCheck();

  void eventUpdate(const CIntEvent&, const C2ParticleData&);

  void eventUpdate(const CGlobEvent&, const CNParticleData&) {}
  
  void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&) {}

  COutputPlugin *Clone() const { return new COPCollDistCheck(*this); }

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

private:

  const COPCollMatrix* ptrCM;
  std::map<std::pair<size_t, EEventType>, C1DHistogram> distList;
  
};

#endif
