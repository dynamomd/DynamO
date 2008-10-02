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

#ifndef CCoordinator_H
#define CCoordinator_H

#include <boost/program_options.hpp>
#include <vector>
#include "../datatypes/pluginpointer.hpp"
#include "engine/engine.hpp"

class CCoordinator
{
private:
public:
  CCoordinator():
    Engine(NULL)
  {}

  boost::program_options::variables_map& parseOptions(int, char*[]);

  void initialise();

  void runSimulation();

  void outputData();

  void outputConfigs();

  void signal_handler(int);

private:
  boost::program_options::variables_map vm;
  //boost::scoped_array<CSimulation> Simulations;
  //bool peekmode;

  smrtPlugPtr<CEngine> Engine;
};

#endif
