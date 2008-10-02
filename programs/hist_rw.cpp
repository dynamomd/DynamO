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
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <boost/program_options.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include "../src/extcode/xmlParser.h"
#include "../src/base/is_exception.hpp"
#include <fenv.h>

using namespace std;
using namespace boost;

static const size_t NGamma = 1;
//Set in the main function
static long double alpha;
static long double minErr;
static size_t dataWindow;
static size_t NStepsPerStep = 0;
static boost::program_options::variables_map vm;

struct SimData;

std::vector<SimData> SimulationData;
//Backup the full list of simulations when piecemeal
std::vector<SimData> TrueList;

std::vector<std::pair<long double, long double> > densOStates;

struct SimData
{
  SimData(std::string nfn):fileName(nfn), logZ(0.0), new_logZ(0.0), refZ(false)
  {
    namespace io = boost::iostreams;
    XMLNode xMainNode;
    if (std::string(fileName.end() - 4, fileName.end()) == ".xml")
      {
	if (!boost::filesystem::exists(fileName))
	  I_throw() << "Could not open XML configuration file";
	
	xMainNode=XMLNode::openFileHelper(fileName.c_str(), "OutputData");
      }
    else if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
      {
	if (!boost::filesystem::exists(fileName))
	  I_throw() << "Could not open XML configuration file";
	
	io::filtering_istream inputFile;
	inputFile.push(io::bzip2_decompressor());
	inputFile.push(io::file_source(fileName));
	//Copy file to a string
	std::string line, fileString;
	
	while(getline(inputFile,line)) 
	  {
	    fileString.append(line);
	    fileString.append("\n");
	  }
	
	XMLNode tmpNode = XMLNode::parseString(fileString.c_str());
	xMainNode = tmpNode.getChildNode("OutputData");
      }
    else
      I_throw() << "Unrecognised extension for input file";

    //Now navigate to the histogram and load the data
    std::istringstream HistogramData(xMainNode.getChildNode("EnergyHist")
				     .getChildNode("WeightHistogram").getText());
    
    //Load gamma here
    if (xMainNode.hasChild("Energy"))
      gamma.push_back(-1.0 / boost::lexical_cast<long double>
		      (xMainNode.getChildNode("Energy").getChildNode("T")
		       .getAttribute("val")));
    else
      gamma.push_back(-1.0 / boost::lexical_cast<long double>
		      (xMainNode.getChildNode("KEnergy").getChildNode("T")
		       .getAttribute("val")));

    std::vector<long double> tmpData;

    tmpData.resize(NGamma+1);

    while (HistogramData >> tmpData[0])
      {
	for (size_t i = 1; i < NGamma+1; ++i)
	  HistogramData >> tmpData[i];

	data.push_back(tmpData);
	tmpData.clear();
	tmpData.resize(NGamma+1);	  
      }
    
  }
  
  bool operator<(const SimData& d2) const
  { return gamma[0] < d2.gamma[0]; }

  bool operator>(const SimData& d2) const
  { return gamma[0] > d2.gamma[0]; }

  std::string fileName;
  std::vector<long double> gamma;
  long double logZ;
  long double new_logZ;
  bool refZ;
  std::vector<std::vector<long double> > data;
};

long double
getlogZ(size_t ID)
{
  size_t bottom = 0, top = SimulationData.size();

  const std::vector<long double>& gamma = SimulationData[ID].gamma;

  if (dataWindow)
    {
      int bID = ID - dataWindow;
      if (bID < 0)
	bottom = 0;
      else 
	bottom = static_cast<size_t>(bID);

      top = ID + dataWindow;
      if (top > SimulationData.size())
	top = SimulationData.size();
    }
    
  long double sum1 = 0.0;

  for (size_t i = bottom; i < top; ++i)
    BOOST_FOREACH(const std::vector<long double>& simdat, SimulationData[i].data)
      {
	long double sum2 = 0.0;
	for (size_t j = bottom; j < top; ++j)
	  {
	    long double dot = 0.0;
	    for (size_t i = 0; i < NGamma; ++i)
	      dot += (SimulationData[j].gamma[i] - gamma[i]) * simdat[i];
	    sum2 += exp(dot - SimulationData[j].logZ);
	  }
	
	sum1 += simdat.back() / sum2;
      }
  
  return log(sum1);
}

void
solveWeights()
{
  std::cout << "##################################################\n";
  std::cout << "Solving for Z's, Error below\n";

  double err = 0.0;

  do
    {
      for (size_t i = NStepsPerStep; i != 0; --i)
	{
	  for (size_t i = 0; i < SimulationData.size(); ++i)
	    if (!SimulationData[i].refZ)
	      SimulationData[i].new_logZ = getlogZ(i);
	  
	  BOOST_FOREACH(SimData& dat, SimulationData)
	    dat.logZ = dat.new_logZ;
	}

      //Now the error checking run
      err = 0.0;
      for (size_t i = 0; i < SimulationData.size(); ++i)
	if (!SimulationData[i].refZ)
	  {
	    SimulationData[i].new_logZ = getlogZ(i);
	    long double tmperr = fabs ((SimulationData[i].new_logZ - SimulationData[i].logZ) / SimulationData[i].new_logZ);
	    if (tmperr > err)
	      err = tmperr;
	  }
      
      //Now copy over the new values
      BOOST_FOREACH(SimData& dat, SimulationData)
	dat.logZ = dat.new_logZ;      

      printf("\r%E", err);
      fflush(stdout);
    }
  while(err > minErr);

  std::cout << "\nIteration complete\n";
}


void calcDensityOfStates()
{
  densOStates.clear();

  typedef std::pair<long double, long double> localpair;

  std::map<long double, long double> accumilator;
  std::cout << "##################################################\n";
  std::cout << "Density of states\n";
  
  //Box up the internal enery histograms
  BOOST_FOREACH(const SimData& dat, SimulationData)
    BOOST_FOREACH(const std::vector<long double>& simdat, dat.data)
    accumilator[simdat.front()] += simdat.back();

  {
    long double sum = 0.0;
    BOOST_FOREACH(const localpair& dat, accumilator)
      sum += dat.second;
    
    std::cout << "Total weight of all data = " << sum << "\n";
  }
  
  BOOST_FOREACH(const localpair& dat, accumilator)
    {
      long double sum = 0.0;
      BOOST_FOREACH(const SimData& dat2, SimulationData)
	sum += exp(dat2.gamma[0] * dat.first - dat2.logZ);
      
      densOStates.push_back(std::make_pair(dat.first, dat.second / sum));
    }
}

void outputDensityOfStates()
{
  typedef std::pair<long double, long double> localpair;
  
  std::fstream of("StateDensity.out",  ios::trunc | ios::out);

  of << std::setprecision(std::numeric_limits<long double>::digits10);

  BOOST_FOREACH(const localpair& dat, densOStates)
    of << dat.first << " " << dat.second << "\n";

  of.close();
}

void outputLogZ()
{
  typedef std::pair<long double, long double> localpair;
  
  std::fstream of("logZ.out",  ios::trunc | ios::out);

  of << std::setprecision(std::numeric_limits<long double>::digits10);

  BOOST_FOREACH(const SimData& dat, SimulationData)
    of << dat.gamma[0] << " " << dat.logZ << "\n";

  of.close();
}

void outputMoments()
{
  std::cout << "##################################################\n";
  std::cout << "Calculating  moments\n";
  typedef std::pair<long double, long double> localpair;
  
  std::vector<long double> betas;

  BOOST_FOREACH(const SimData& dat, SimulationData)
    betas.push_back(dat.gamma[0]);

  long double step = (betas.front() - betas.back()) / 1000;

  std::vector<std::pair<double,double> > Cv;

  {
    std::fstream Eof("Energy.out",  ios::trunc | ios::out);
    std::fstream E2of("Energy2.out",  ios::trunc | ios::out);
    std::fstream Cvof("Cv.out",  ios::trunc | ios::out);
    
    Eof << std::setprecision(std::numeric_limits<long double>::digits10);
    E2of << std::setprecision(std::numeric_limits<long double>::digits10);
    Cvof << std::setprecision(std::numeric_limits<long double>::digits10);
    
    
    
    for (long double Beta = betas.back(); Beta >= betas.front(); Beta += step)
      {
	//Calc Z
	long double Z = 0.0;
	BOOST_FOREACH(const localpair& dat, densOStates)
	  Z += exp (log(dat.second) + Beta * dat.first);
	
	Z = log(Z);
	
	//Now calc the normalisation, or first moment
	long double Norm = 0.0;
	long double Eavg = 0.0;
	long double E2avg = 0.0;
	
	BOOST_FOREACH(const localpair& dat, densOStates)
	  {
	    long double temp = exp(log(dat.second) + Beta * dat.first - Z);
	    Norm += temp;
	    Eavg += temp * dat.first;
	    E2avg += temp * dat.first * dat.first;
	  }
	
	Eavg /= Norm;
	E2avg /= Norm;
	
	Eof << -1.0/Beta << " " << Eavg << "\n";
	E2of << -1.0/Beta << " " << E2avg  << "\n";
	Cv.push_back(std::make_pair(-1.0/Beta, Beta * Beta * (E2avg - Eavg * Eavg)));
	Cvof << -1.0/Beta << " " << Beta * Beta * (E2avg - Eavg * Eavg)  << "\n";
      }
    Eof.close();
    E2of.close();
    Cvof.close();
  }

  {
    std::fstream CvMax("Cvmax.out",  ios::trunc | ios::out);
    std::fstream CvMin("Cvmin.out",  ios::trunc | ios::out);
    
    CvMax << std::setprecision(std::numeric_limits<long double>::digits10);
    CvMin << std::setprecision(std::numeric_limits<long double>::digits10);
    
    //Now analyse the Cv for its turning points etc
    double olddt = (Cv.begin()+1)->second - Cv.begin()->second;
    for (std::vector<std::pair<double,double> >::const_iterator iPtr = Cv.begin()+2; iPtr != Cv.end(); ++iPtr)
      {
	double thisdt = iPtr->second - (iPtr - 1)->second;
	if (std::signbit(olddt) != std::signbit(thisdt))
	  if (std::signbit(olddt))
	    CvMin << (iPtr-1)->first << " " << (iPtr-1)->second << "\n";
	  else
	    CvMax << (iPtr-1)->first << " " << (iPtr-1)->second << "\n";
	olddt = thisdt;
      }
    
    CvMax.close();
    CvMin.close();
  }
  
}

void piecemealRun(size_t lower_ID, size_t upper_ID, size_t refZ)
{
  SimulationData.clear();

  for (size_t i = lower_ID; i < upper_ID; ++i)
    {
      SimulationData.push_back(TrueList[i]);
      //Correctly set what is the reference Z
      SimulationData.back().refZ = (i == refZ);
    }

  solveWeights();

  //Copy the data back now
  for (size_t i = lower_ID; i < upper_ID; ++i)
    {
      TrueList[i] = SimulationData[i - lower_ID];
      TrueList[i].refZ = false;
    }

  SimulationData.clear();
}

void solveWeightsPiecemeal()
{
  //Backup the full list of simulations
  TrueList = SimulationData;
  SimulationData.clear();
	
  size_t logZID = 0;
	
  //Find the zero or reference logZ
  for (; logZID < TrueList.size(); ++logZID)
    if (TrueList[logZID].refZ)
      break;
	
  //Start moving up first
  int lower_ID = logZID, 
    upper_ID = lower_ID + vm["piecemeal"].as<size_t>();
		
  while (upper_ID <= TrueList.size())
    {
      std::cout << "\nPiece " << lower_ID << " to " << upper_ID - 1;
      piecemealRun(lower_ID, upper_ID, lower_ID);
      //Move again
      lower_ID = upper_ID - 1;
      upper_ID = lower_ID + vm["piecemeal"].as<size_t>();
    }

  if (lower_ID != TrueList.size())
    {
      //One last one
      upper_ID = TrueList.size();
      std::cout << "\nLast Piece " << lower_ID << " to " << upper_ID - 1;
      piecemealRun(lower_ID, upper_ID, lower_ID);	    
    }

  //Start moving down
  upper_ID = logZID + 1, 
    lower_ID = upper_ID - vm["piecemeal"].as<size_t>();

  while (lower_ID >= 0)
    {
      std::cout << "\nPiece " << lower_ID << " to " << upper_ID - 1;
      piecemealRun(lower_ID, upper_ID, upper_ID - 1);
      //Move again
      upper_ID = lower_ID + 1;
      lower_ID = upper_ID - vm["piecemeal"].as<size_t>();
    }

  if (upper_ID != 0)
    {
      //One last one
      lower_ID = 0;
      std::cout << "\nLast Piece " << lower_ID << " to " << upper_ID - 1;
      piecemealRun(lower_ID, upper_ID, upper_ID - 1);
    }
	
  //Restore the SimulationData list
  SimulationData = TrueList;
}

void solveWeightsRollerPiecemeal()
{
  std::cout << "##################################################\n";
  std::cout << "Piecemeal roller solution";

  BOOST_FOREACH(SimData& dat, SimulationData)
    dat.refZ = false;

  //Because of this we must roll forward
  SimulationData.front().refZ = true;
  SimulationData.front().new_logZ = SimulationData.front().logZ;

  if (vm["piecemeal"].as<size_t>() + 1 > SimulationData.size())
    I_throw() << "Not enough data for that piecemeal roller solution";

  std::list<SimData> StartList;
  std::list<SimData> FinalList;
  
  for (size_t i = 0; i < SimulationData.size(); ++i)
    StartList.push_back(SimulationData[i]);

  SimulationData.clear();

  for (size_t i = 0; i < vm["piecemeal"].as<size_t>(); ++i)
    {
      SimulationData.push_back(StartList.front());
      StartList.pop_front();
    }

  std::cout <<"Initial run takes the longest\n";

  long double oldval = SimulationData.back().logZ;
  
  solveWeights();

  long double diff =  SimulationData.back().logZ - oldval;

  while (!StartList.empty())
    {
      std::cout <<"################ COUNTDOWN " << StartList.size()  << "\n";
      
      for (size_t i = 0; i < SimulationData.size()/2; ++i)
	{
	  SimulationData[i].refZ = true;
	  SimulationData[i].new_logZ = SimulationData[i].logZ;
	}
      
      FinalList.push_back(SimulationData.front());

      for (size_t i = 0; i < SimulationData.size()-1; ++i)
	  SimulationData[i] = SimulationData[i+1];
      
      SimulationData.back() = StartList.front();
      StartList.pop_front();

      oldval = SimulationData.back().logZ;
      SimulationData.back().logZ += diff;

      solveWeights();
      diff = SimulationData.back().logZ - oldval;
    }

  for (size_t i = 0; i < SimulationData.size(); ++i)
    FinalList.push_back(SimulationData[i]);
  
  SimulationData.clear();
  
  BOOST_FOREACH(const SimData& dat, FinalList)
    SimulationData.push_back(dat);
}

int
main(int argc, char *argv[])
{
  std::cout << "dynahist_rw  Copyright (C) 2008  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n\n";

  //This is so the program crashes out when floating point errors occur
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  try {
    namespace po = boost::program_options;
    
    boost::program_options::options_description systemopts("Program Options");
    
    systemopts.add_options()
      ("help", "Produces this message")   
      ("data-file", po::value<std::vector<std::string> >(), "Specify a config file to load, or just list them on the command line")
      ("alpha", po::value<long double>()->default_value(1), "A fraction of the difference between the old and new logZ's to use, use to stop divergence")
      ("NSteps,N", po::value<size_t>()->default_value(10), "Number of steps to take before testing the error and spitting out the current vals")
      ("load-logZ", po::value<std::string>(), "Loads the logZ's from a previous run, note! It does this by ordering the temperatures and adding in order, do not change anything you do!")
      ("min-err", po::value<long double>()->default_value(1.0e-5), "The minium error allowed before the loop terminates")
      ("data-window,h", po::value<size_t>()->default_value(0), "The +- simulations around a simulation data to use in calculating Z (0 is all)")
      ("piecemeal,p", po::value<size_t>(), "If set then the solution is split into chunks of the specified size and solved piecemeal before being combined for the global solution, used to get a fast initial guess which can be very near the final solution")
      ("fast-piecemeal", "Better piecemeal derivative, requires piecemeal be set")
      ;

    boost::program_options::positional_options_description p;
    p.add("data-file", -1);
    
    boost::program_options::store(po::command_line_parser(argc, argv).
				  options(systemopts).positional(p).run(), vm);
    boost::program_options::notify(vm);
    
    if (vm.count("help") || !vm.count("data-file")) 
      {
	I_throw() << "Usage : dynahist_rw <OPTION>...<data-file(s)>\n"
		  << "Determines the weighting functions for the histograms\n"
		  << systemopts << "\n";
      }

    alpha = vm["alpha"].as<long double>();
    NStepsPerStep = vm["NSteps"].as<size_t>();
    minErr = vm["min-err"].as<long double>();
    dataWindow = vm["data-window"].as<size_t>();

    //Data load
    BOOST_FOREACH(std::string fileName, vm["data-file"].as<std::vector<std::string> >())
      SimulationData.push_back(SimData(fileName));
    
    //Sort by temperature
    std::sort(SimulationData.begin(), SimulationData.end());

    //Output ordered sims
    std::cout << "##################################################\n";
    BOOST_FOREACH(const SimData& dat, SimulationData)
      std::cout << dat.fileName << " NData = " << dat.data.size() << " gamma[0] = " << dat.gamma[0] << "\n";

    if (vm.count("load-logZ"))
      {
	std::fstream logZin(vm["load-logZ"].as<std::string>().c_str(), std::ios::in);

	long double tmp;

	BOOST_FOREACH(SimData& dat, SimulationData)
	  {
	    logZin >> tmp;
	    logZin >> tmp;
	    dat.logZ = tmp;
	    dat.new_logZ = tmp;
	  }
	logZin.close();

	bool done = false;
	BOOST_FOREACH(SimData& dat, SimulationData)
	  if (dat.logZ == 0.0)
	    {
	      dat.refZ = true;
	      done = true;
	      break;
	    }

	if (!done)
	  //set the reference system about midway to reduce over/underflows
	  (SimulationData.begin() + (SimulationData.size()/2) )->refZ = true;	
      }
    else
      //set the reference system about midway to reduce over/underflows
      (SimulationData.begin() + (SimulationData.size()/2) )->refZ = true;

    if (vm.count("piecemeal"))
      {
	if (vm.count("fast-piecemeal"))
	  solveWeightsRollerPiecemeal();
	else 
	  solveWeightsPiecemeal();
      }
    else
      solveWeights();
    
    std::cout << "##################################################\n";
    BOOST_FOREACH(const SimData& dat, SimulationData)
      std::cout << dat.fileName << " logZ = " << dat.logZ << "\n";
        
    //Now start the fun
    outputLogZ();
    calcDensityOfStates();
    outputDensityOfStates();
    outputMoments();
    
  }
  catch (DYNAMO::Exception& cep)
    {
      fflush(stdout);
      std::cerr << cep.what() << "\nMAIN: Reached Main Error Loop\n";
      return 1;
    }
}
