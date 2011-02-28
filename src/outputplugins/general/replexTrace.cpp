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

#include "replexTrace.hpp"
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include "../../dynamics/include.hpp"
#include "../../extcode/xmlwriter.hpp"

OPReplexTrace::OPReplexTrace(const DYNAMO::SimData* t1, const XMLNode&):
  OutputPlugin(t1,"ReplexTrace"),
  filename(std::string("ReplexTrace.tmp.") 
	   + boost::lexical_cast<std::string>(rand()))
{
  while (boost::filesystem::exists(filename))
    {
      filename = std::string("ReplexTrace.tmp.") 
	+ boost::lexical_cast<std::string>(rand());
    }
  
  tmpfile.open(filename.c_str(), std::ios::in | std::ios::out 
	       | std::ios::trunc);

  if (!tmpfile.is_open())
    M_throw() << "Could not open temporary file!";
}

OPReplexTrace::~OPReplexTrace()
{
  //Clean up temporary files
  tmpfile.close();
  boost::filesystem::remove(filename);
}

OPReplexTrace::OPReplexTrace(const OPReplexTrace& cop2):
  OutputPlugin(cop2),
  filename(std::string("./ReplexTrace.tmp.") 
	   + boost::lexical_cast<std::string>(rand()))
{
  while (boost::filesystem::exists(filename))
    {
      filename = std::string("./ReplexTrace.tmp.") 
	+ boost::lexical_cast<std::string>(rand());
    }

  tmpfile.open(filename.c_str(), std::ios::in | std::ios::out 
	       | std::ios::trunc);

  if (!tmpfile.is_open())
    M_throw() << "Could not open temporary file!";

  //Copy the file stream
  cop2.tmpfile.seekg (0, std::ios::beg);
  std::copy(std::istreambuf_iterator<char>(cop2.tmpfile),
	    std::istreambuf_iterator<char>(),
	    std::ostreambuf_iterator<char>(tmpfile));
}

void 
OPReplexTrace::initialise() 
{ 
  if (!(tmpfile.is_open()))
    M_throw() << "OPReplexTrace temp file unopened!";
}

void 
OPReplexTrace::changeSystem(OutputPlugin* OPP)
{
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<OPReplexTrace*>(OPP) == NULL)
    M_throw() << "Not the correct plugin to change System with";
#endif

  addPoint();
  static_cast<OPReplexTrace*>(OPP)->addPoint();

  std::swap(Sim, static_cast<OPReplexTrace*>(OPP)->Sim);

  addPoint(); 
  static_cast<OPReplexTrace*>(OPP)->addPoint();  
}

void 
OPReplexTrace::addPoint()
{
  const boost::array<double,3>& ensembleVals(Sim->ensemble->getReducedEnsembleVals());
    
  tmpfile << Sim->dSysTime / Sim->dynamics.units().unitTime() 
	  << " " 
	  << ensembleVals[0] << "," << ensembleVals[1] << "," << ensembleVals[2] << "\n";
}

void 
OPReplexTrace::output(xml::XmlStream& XML)
{
  addPoint();

  XML << xml::tag("ReplexTrace")
      << xml::chardata();
  
  tmpfile.seekg (0, std::ios::beg);
  std::copy(std::istreambuf_iterator<char>(tmpfile),
	    std::istreambuf_iterator<char>(),
	    std::ostreambuf_iterator<char>(XML.getUnderlyingStream()));  

  XML << xml::endtag("ReplexTrace");
}
