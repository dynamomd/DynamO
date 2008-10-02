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

#ifndef COPTRAJECTORY_H
#define COPTRAJECTORY_H

#include "outputplugin.hpp"
#include <string>
#include <list>

//Must use a collHistory list to allow copying and saving of state,
//Cannot write directly to file!
class mempiece;

class COPTrajectory : public COutputPlugin
{  
 public:
  COPTrajectory(DYNAMO::SimData*);
  ~COPTrajectory();

  void output(xmlw::XmlStream &); 

  void collisionUpdate(const CIntEvent &, const CIntEventData &);

  virtual COutputPlugin *Clone() const
    { return new COPTrajectory(*this); }

  void setFilename(std::string fn)
    { fileName = fn; }
  
  void setSimulationType(unsigned int a)
    { simType = a; }

  void setTotalCollCount(unsigned long);

 private:

  static std::list<mempiece> collHistory;
  std::list<mempiece>::iterator currentPos;
  std::string fileName;
  unsigned int simType;
};

#endif
