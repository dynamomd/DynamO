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

#ifndef CEngine_H
#define CEngine_H

#include <boost/program_options.hpp>
#include <boost/scoped_array.hpp>
#include "../../simulation/simulation.hpp"
#include "../../extcode/threadpool.hpp"

class CIPCompression;

class CEngine
{
public:
  CEngine(const boost::program_options::variables_map&,
	  std::string, std::string);
  
  virtual ~CEngine() {}

  virtual void initialisation() = 0;

  virtual void finaliseRun() = 0;

  virtual void forceShutdown() = 0;

  virtual void printStatus() = 0;

  virtual void runSimulation() = 0;

  virtual void outputData() = 0;

  virtual void outputConfigs() = 0;
  
  virtual void peekData() = 0;
  
  static void getCommonOptions(boost::program_options::options_description&);
  
protected:
  virtual void preSimInit();

  virtual void setupSim(CSimulation &, const std::string);

  virtual void postSimInit(CSimulation&);

  const boost::program_options::variables_map& vm;
  //std::vector<CIPCompression*> CompPlugs;
  CThreadPool threads;
  std::string configFormat;
  std::string outputFormat;
};

#endif
