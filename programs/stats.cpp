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
#include <string>
#include <math.h>
#include <iostream>
#include <boost/program_options.hpp>
#include <iomanip>

using namespace boost;
using namespace std;
namespace po = boost::program_options;

#include "stats.hpp"
#include "../base/is_exception.hpp"

int main(int ac, char* av[])
{
  try 
    {      
      po::options_description desc("Allowed options");
      desc.add_options()
	("help", "produce help message")
	("input-file", po::value< std::vector<std::string > >(), "input files for processing")
	("tag,t", po::value<std::vector<std::string> >(), "XML tag name bounding the data to be processed")
	("attr", po::value<std::string>()->default_value("val"), "Attribute to average")
	("cut-mode,c", "Just cut the column out")
	("graph-mode", "process columns of data enclosed in xml tags, last column will be averaged and stdev calculated")
	("info-mode,i", "Gives information on the output file")
	;
      
      po::positional_options_description p;
      p.add("input-file", -1);
      
      po::variables_map vm;
      po::store(po::command_line_parser(ac, av).
		options(desc).positional(p).run(), vm);
      po::notify(vm);
      
      if (vm.count("help")) {
	std::cerr << "Usage: options_description [options]\n";
	std::cerr << desc;
	return 0;
      }

      std::cout << std::setprecision(std::numeric_limits<Iflt>::digits10);

      vector<string> files,tags;

      if (vm.count("input-file"))
        files = vm["input-file"].as<std::vector<std::string> >();
      else
	  I_throw() << "No input files specified...exiting";

      if (vm.count("tag"))
        tags = vm["tag"].as<std::vector<std::string> >();
      else
	  I_throw() << "No tags specified...exiting";

      XMLNode xMainNode, xBrowseNode;
      vector< vector<Iflt> > graphdata;
      vector< Iflt > listdata;

      //next we have a list of files to process
      for (vector<string>::const_iterator iPtr = files.begin(); iPtr != files.end(); iPtr++)
	{
	  
	  //std::cerr << "Processing file " << *iPtr << std::endl;
	  io::filtering_istream coutputFile;
	  coutputFile.push(io::bzip2_decompressor());
	  coutputFile.push(io::file_source(*iPtr));

	  if (!(coutputFile.component<1,io::file_source>()->is_open()))
	    I_throw() << "Could not open data file!";
	  
	  string line,fileString;
	  while(getline(coutputFile,line)) 
	    {
	      fileString.append(line);
	      fileString.append("\n");
	    }
	  
	  //std::cerr << "Parsing " << std::endl;
	  xMainNode = XMLNode::parseString(fileString.c_str());
	  xBrowseNode = xMainNode;
	  	  
	  XMLNode xSubNode;

	  //std::cerr << "Navigating " << std::endl;
	  //Navigate to the correct tag
	  try {
	    for (vector<string>::const_iterator jPtr = tags.begin(); jPtr != tags.end(); jPtr++)
	      {
		xSubNode = xBrowseNode.getChildNode(jPtr->c_str());
		xBrowseNode = xSubNode;
	      }
	    
	    //Ready to average!
	    if (vm.count("graph-mode"))
	      {
		//Make the graphdata big enough
		if (graphdata.size() == 0)
		  graphdata = vector<vector<Iflt> >(xBrowseNode.nText());
		
		//Push on each line
		for (int i = 0; i < xBrowseNode.nText(); i++)
		  (graphdata[i]).push_back(atof(xBrowseNode.getText(i)));
		
	      } 
	    else if (vm.count("cut-mode"))
	      {
		for (int i = 0; i < xBrowseNode.nText(); i++)
		  std::cout << xBrowseNode.getText(); 
	      }
	    else
	      listdata.push_back(atof(xBrowseNode.getAttribute(vm["attr"].as<std::string>().c_str())));
	  } 
	  catch (...)
	    {
	      std::cerr << "Failed in navigation of the file " << *iPtr << endl;
	    }
	}

      //Ok now we have them as STL vectors of data, lets work
      if (vm.count("cut-mode"))
	{}
      else if (vm.count("graph-mode"))
	{} 
      else
	{
	  Iflt avg = 0.0;
	  for (vector<Iflt>::const_iterator iPtr = listdata.begin();
	       iPtr != listdata.end(); iPtr++)
	    avg += *iPtr;
	  avg /= listdata.size();

	  Iflt SD = 0.0;
	  for (vector<Iflt>::const_iterator iPtr = listdata.begin();
	       iPtr != listdata.end(); iPtr++)
	    SD += ((*iPtr) - avg)*((*iPtr) - avg);
	  
	  SD /= listdata.size ();
	  
	  SD = sqrt(SD);
	  
	  cout << "Samples " << listdata.size() << "\n";
	  cout << "Average " << avg << "\n";
	  cout << "SD      " << SD << "\n";
	}
      
    }
  catch(exception& e)
    {
      cerr << e.what() << "\n";
      return 1;
    }    

  return 0;
}
