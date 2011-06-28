/*  dynamo:- Event driven molecular dynamics simulator 
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

#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <cmath>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <fenv.h>
#include <buildinfo.hpp>

using namespace std;

struct data
{
  double T;
  double upSims;
  double downSims;
  double f;
};

int main(int argc, char* argv[])
{
  //The following macro converts the GITHASH define to a C style
  //string, the boost build system won't let us define strings on the
  //command line.
#define VALUE_TO_STRING(val) #val
#define STR(val) VALUE_TO_STRING(val)
  std::cout << "dynareplex_opt Copyright (C) 2011  Marcus N Campbell Bannerman\n"
            << "This program comes with ABSOLUTELY NO WARRANTY.\n"
            << "This is free software, and you are welcome to redistribute it\n"
            << "under certain conditions. See the licence you obtained with\n"
            << "the code\n"
               "Git Checkout Hash " << STR(GITHASH) << "\n\n";
  //This is so the program crashes out when floating point errors occur
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

    namespace po = boost::program_options;
    
    boost::program_options::options_description systemopts("Program Options");
    
    systemopts.add_options()
      ("help", "Produces this message")   
      ("three-point,3", "Use 3 point method for the derivatives, can stop floating point errors in bad statistics runs")
      ("data-file", po::value<std::string>()->default_value("replex.dat"), "File containing 5 columns, 1st is temperature, 4th and 5th are the low-high and high-low counts")
      ("alpha", po::value<double>()->default_value(1.0), "A fraction of the difference between the old and new T's to use. debugging use really")
      ("NSims,S", po::value<size_t>(), "Output optimised temperatures for a different number of sims, default is the number of input sims")
      ("configmod-commands,c", "For output print the commands to is_configmod to set the temperatures")
      ;

    boost::program_options::variables_map vm;
    boost::program_options::store(po::command_line_parser(argc, argv).
				  options(systemopts).run(), vm);
    boost::program_options::notify(vm);
    
    if (vm.count("help") || !vm.count("data-file")) 
      {
	std::cout << "Usage : dynareplex_opt <OPTIONS>\n"
		  << "Gives the next iteration of best temperatures for replica exchange\n"
		  << systemopts << "\n";
	exit(1);
      }

  std::vector<data> mydata;
 
  {
    data tmpdata;
    
    fstream inputf(vm["data-file"].as<std::string>().c_str(), ios::in);
    
    while (inputf >> tmpdata.T >> tmpdata.upSims >> tmpdata.upSims >> tmpdata.upSims >> tmpdata.downSims)
      mydata.push_back(tmpdata);
    
    inputf.close();
  }

  BOOST_FOREACH(data& dat, mydata)
    dat.f = dat.downSims / (dat.upSims + dat.downSims);

  {
    fstream ff("f.out", ios::out|ios::trunc);
    
    for (size_t i = 0; i < mydata.size(); ++i)
      {      
	ff << mydata[i].T << " " << mydata[i].f << "\n";
      }
    ff.close();
  }

  size_t N = mydata.size();

  std::vector<double> dT;

  for (size_t i = 1; i < N; ++i)
    dT.push_back(mydata[i].T - mydata[i-1].T);

  std::vector<double> dfdT;


  for (size_t i = 1; i < N; ++i)
    if (((N-i) > 1) && vm.count("three-point"))
      //Central 3 point method
      dfdT.push_back((mydata[i+1].f - mydata[i-1].f)
		     / (mydata[i+1].T - mydata[i-1].T));
    else
      //First order fall back (two data points or not enabling anything)
      dfdT.push_back((mydata[i].f - mydata[i-1].f) / dT[i-1]);

  {
    fstream dfdtf("dfdt.out", ios::out|ios::trunc);
    for (size_t i = 0; i < dT.size(); ++i)
      dfdtf << mydata[i].T << " " << dfdT[i] << "\n"
	    << mydata[i+1].T << " " << dfdT[i] << "\n";

    dfdtf.close();
  }


  double norm = 0.0;
  for (size_t i = 0; i < dT.size(); ++i)
    norm += sqrt(dfdT[i] * dT[i]);
  
  norm = 1.0 / norm;

  {
    fstream etaf("eta.out", ios::out|ios::trunc);
    for (size_t i = 0; i < dT.size(); ++i)
      {      
	double val = sqrt(dfdT[i]/ dT[i]);
	etaf << mydata[i].T << " " << norm * val << "\n"
	     << mydata[i+1].T << " " << norm * val << "\n";
      }

    etaf.close();
  }

  size_t index = 0;

  double sumdiff = 0.0;

  size_t NSims = mydata.size();
  if (vm.count("NSims"))
    NSims = vm["NSims"].as<size_t>();

  for (size_t i = 0; i < NSims - 1; ++i)
    {
      double next_ival = norm * sqrt(dfdT[index] * dT[index]);

      double target_prob = static_cast<double>(i) / NSims;

      //keep moving forward
      while (sumdiff + next_ival < target_prob)
	{
	  sumdiff += next_ival;
	  ++index;
	  next_ival = norm * sqrt(dfdT[index] * dT[index]);
	}
      
      double T = ((target_prob - sumdiff) * dT[index] / next_ival) + mydata.at(index).T;

      double dT = vm["alpha"].as<double>() * (T - mydata.at(i).T);

      if (!vm.count("configmod-commands"))
	std::cout << i << " " 
		  << mydata.at(i).T + dT 
		  << " " << dT
		  << "\n";
      else
	std::cout << "is_configmod -T " << mydata.at(i).T + dT << " config." 
		  << i << ".end.xml.bz2 -o config." << i << ".end.xml.bz2\n";
    }
  
      if (!vm.count("configmod-commands"))
	std::cout << NSims - 1 << " " << mydata.back().T
		  << " 0\n";
      else
	std::cout << "is_configmod -T " << mydata.back().T << " config." 
		  << NSims - 1 << ".end.xml.bz2 -o config." << NSims - 1 
		  << ".end.xml.bz2\n";       
}
