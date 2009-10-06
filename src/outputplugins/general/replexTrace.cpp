/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

COPReplexTrace::COPReplexTrace(const DYNAMO::SimData* t1, const XMLNode&):
  COutputPlugin(t1,"ReplexTrace"),
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
    D_throw() << "Could not open temporary file!";
}

COPReplexTrace::~COPReplexTrace()
{
  //Clean up temporary files
  tmpfile.close();
  boost::filesystem::remove(filename);
}

COPReplexTrace::COPReplexTrace(const COPReplexTrace& cop2):
  COutputPlugin(cop2),
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
    D_throw() << "Could not open temporary file!";

  //Copy the file stream
  cop2.tmpfile.seekg (0, std::ios::beg);
  std::copy(std::istreambuf_iterator<char>(cop2.tmpfile),
	    std::istreambuf_iterator<char>(),
	    std::ostreambuf_iterator<char>(tmpfile));
}

void 
COPReplexTrace::initialise() 
{ 
  if (!(tmpfile.is_open()))
    D_throw() << "COPReplexTrace temp file unopened!";
}

void 
COPReplexTrace::changeSystem(COutputPlugin* OPP)
{
#ifdef DYNAMO_DEBUG
  if (dynamic_cast<COPReplexTrace*>(OPP) == NULL)
    D_throw() << "Not the correct plugin to change System with";
#endif

  addPoint();
  static_cast<COPReplexTrace*>(OPP)->addPoint();

  std::swap(Sim, static_cast<COPReplexTrace*>(OPP)->Sim);

  addPoint(); 
  static_cast<COPReplexTrace*>(OPP)->addPoint();  
}

void 
COPReplexTrace::addPoint()
{
  tmpfile << Sim->dSysTime << " " << Sim->simID << "\n";
}

void 
COPReplexTrace::output(xmlw::XmlStream& XML)
{
  addPoint();

  XML << xmlw::tag("ReplexTrace")
      << xmlw::chardata();
  
  tmpfile.seekg (0, std::ios::beg);
  std::copy(std::istreambuf_iterator<char>(tmpfile),
	    std::istreambuf_iterator<char>(),
	    std::ostreambuf_iterator<char>(XML.getUnderlyingStream()));  

  XML << xmlw::endtag("ReplexTrace");
}
