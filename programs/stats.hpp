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
#ifndef data_obj_H
#define data_obj_H

#include <string>
#include <ostream>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include "../base/constants.hpp"
#include "../extcode/xmlParser.h"
#include <map>
#include <vector>
#include <algorithm>

using namespace boost;
namespace io = boost::iostreams;

class data_object
{
 public:

  void addFiles(std::vector<std::string> files)
    {
      for (std::vector<std::string>::iterator iPtr = files.begin();
	   iPtr != files.end(); iPtr++)
	addFile(*iPtr);
    }
  
  void addFile(std::string fileName)
    {       
      io::filtering_istream coutputFile;
      coutputFile.push(io::bzip2_decompressor());
      coutputFile.push(io::file_source(fileName));
       
      { //Scope deletes the strings when done
	
	//Copy file to a string
	std::string line, fileString;
	while(getline(coutputFile,line)) 
	  {
	    fileString.append(line);
	    fileString.append("\n");
	  }
	
	//Parse the string
	xmlData[fileName] = fileString;	
      }      
    }

  void recursiveProcess(std::string fileName, XMLNode & XML)
    {
      std::cerr << "recursive processing in tag " << XML.getName() 
		<< " with " << XML.nChildNode() << " child nodes\n";
      if (XML.isAttributeSet("val"))
	tags2files[std::string(XML.getName())].push_back(std::string(fileName)); 

      std::cout << XML.getName() << "\n";

      for (int i=0; i<XML.nChildNode(); i++)
	{
	  XMLNode subNode = XML.getChildNode(i);
	  recursiveProcess(fileName, subNode);
	}
    }

  void parseData()
    {
      tags2files.clear();
     
      std::cout << "Parsing data\n";

      for (std::map<std::string, std::string>::iterator iPtr = xmlData.begin();
	   iPtr != xmlData.end(); iPtr++)
	{
	  XMLNode XML = XMLNode::parseString(iPtr->first.c_str());	
	  recursiveProcess(iPtr->first, XML);
	}
    }

  std::vector<std::pair<std::string, long> > tagList() const
    {
      std::vector<std::pair<std::string, long> > tags;
      
      for (std::map<std::string, std::vector<std::string> >::const_iterator iPtr = tags2files.begin();
	   iPtr != tags2files.end(); iPtr++)
	{
	  tags.push_back(std::pair<std::string, long>(iPtr->first, iPtr->second.size()));
	}

      std::sort(tags.begin(), tags.end());
      
      return tags;
    }
  
  std::vector<std::string> tagFiles(std::string tagname) const
    {
      if (tags2files.find(tagname) != tags2files.end())
	return ((tags2files.find(tagname))->second);
      else
	return std::vector<std::string>();
    }

  Iflt recursiveSearch(std::string tag, XMLNode & XML) const
    {
      std::cout << "\n" << XML.getName();
      if (XML.getName() == tag.c_str())
	if (XML.isAttributeSet("val"))
	  return atof(XML.getAttribute("val"));
	else
	  throw std::runtime_error("Tag recursively found but val not set!");

      else
	for (int i=0; i<XML.nChildNode(); i++)
	  {
	    XMLNode subNode = XML.getChildNode(i);
	    Iflt tmp = recursiveSearch(tag, subNode);
	    if (tmp != HUGE_VAL)
	      return tmp;
	  }
      
      return HUGE_VAL;
    }

  Iflt getVal(std::string file, std::string tag) const
    {
      if (xmlData.find(file) == xmlData.end())
	throw std::runtime_error("Couldn't find the file when getting the tag val");

      XMLNode XML = XMLNode::parseString((xmlData.find(file))->second.c_str());

      Iflt tmp = recursiveSearch(tag,XML);
      if (tmp == HUGE_VAL)
	throw std::runtime_error("Couldn't find the tag in the file");

      return tmp;
    }

  std::pair<Iflt,Iflt> avgSDTag(std::string tag, std::string xtag = "", Iflt xval = 0, Iflt tolerance = 0.05) const
    {
      if (tags2files.find(tag) == tags2files.end())
	throw std::runtime_error("Couldn't find the files for the tag");
      
      std::vector<Iflt> vals;
      
      std::vector<std::string> files = tags2files.find(tag)->second;
      
      if (xtag != "")
	//start filtering
	for (std::vector<std::string>::iterator iPtr = files.begin();
	     iPtr != files.end();)
	  {
	    Iflt val = getVal(*iPtr, tag);
	    if (!(fabs((val-xval)/val) < tolerance))
	      files.erase(iPtr);
	    else
	      iPtr++;
	  }
      
      for (std::vector<std::string>::iterator iPtr = files.begin();
	   iPtr != files.end(); iPtr++)
	vals.push_back(getVal(*iPtr, tag));
      
      Iflt avg = 0.0;
      for (std::vector<Iflt>::const_iterator iPtr = vals.begin();
	   iPtr != vals.end(); iPtr++)
	avg += *iPtr;
      avg /= static_cast<Iflt>(vals.size());
      
      Iflt SD = 0.0;
      for (std::vector<Iflt>::const_iterator iPtr = vals.begin();
	   iPtr != vals.end(); iPtr++)
	SD += ((*iPtr) - avg)*((*iPtr) - avg);
      
      SD /= static_cast<Iflt>(vals.size ());
      
      SD = sqrt(SD);

      return std::pair<Iflt,Iflt>(avg,SD);
    }
  
  /*  std::vector<std::pair<Iflt, pair<Iflt, Iflt> > > getData(std::string xval, std::string yval)
    {
      std::vector<std::pair<Iflt, pair<Iflt, Iflt> > > data;
      
      std::map<Iflt,

      }*/

 private:

  std::map<std::string, std::string> xmlData;

  std::map<std::string, std::vector<std::string> > tags2files;

};

#endif
