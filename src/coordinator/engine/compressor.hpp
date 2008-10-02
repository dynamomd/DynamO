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

#ifndef CECompressor_H
#define CECompressor_H

#include "single.hpp"
#include "../../inputplugins/compression.hpp"

class CECompressor: public CESingle
{
public:
  CECompressor(const boost::program_options::variables_map&);

  virtual ~CECompressor() {}

  virtual void finaliseRun();

  static void getOptions(boost::program_options::options_description&);

protected:
  virtual void preSimInit();
  
  virtual void setupSim(CSimulation&, const std::string);

  smrtPlugPtr<CIPCompression> compressPlug;
};

#endif
