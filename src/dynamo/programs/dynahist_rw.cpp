/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <unordered_map>

#include <fenv.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iosfwd>
#include <array>

using namespace std;
using namespace boost;

const size_t NGamma = 1;

//Set in the main function
static long double alpha;
static long double minErr = 1e-16;
static size_t NStepsPerStep = 0;
static boost::program_options::variables_map vm;

long double betaMax;
long double betaMin;

struct SimulationData;

std::vector<SimulationData> SimulationDataData;

struct SimulationData
{
  SimulationData(std::string nfn):fileName(nfn), logZ(0.0), new_logZ(0.0), refZ(false)
  {
    using namespace magnet::xml;

    if (!boost::filesystem::exists(fileName))
      M_throw() << "Could not find the XML file named " << fileName
		<< "\nPlease check the file exists.";

    Document doc(fileName);

    Node mainNode = doc.getNode("OutputData");

    //Find the bin width
    if (!(mainNode.hasNode("EnergyHist")))
      M_throw() << "Could not find the Internal Energy Histogram in output file " << nfn;

    if (!(mainNode.getNode("EnergyHist").hasAttribute("BinWidth")))
      M_throw() << "Could not find the BinWidth attribute in the Internal Energy Histogram";

    if (!(mainNode.getNode("EnergyHist").hasAttribute("T")))
      M_throw() << "Could not find the Temperature attribute in the Internal Energy Histogram";

    binWidth = mainNode.getNode("EnergyHist").getAttribute("BinWidth").as<double>();

    gamma.push_back(-1.0 / (mainNode.getNode("EnergyHist").getAttribute("T").as<long double>()));

    //Load the W factor for each energy
    if (mainNode.getNode("EnergyHist").hasNode("PotentialDeformation"))
      {
	for (magnet::xml::Node node = mainNode.getNode("EnergyHist").getNode("PotentialDeformation").findNode("W"); 
	     node.valid(); ++node)
	  {
	    double energy = node.getAttribute("Energy").as<double>();
	    double Wval = node.getAttribute("Value").as<double>();
	    if (Wval) _W[lrint(energy / binWidth)] = Wval;
	  }
      }

    std::cout << "W for file " << nfn;
    for (std::unordered_map<int, double>::iterator iPtr = _W.begin();
	 iPtr != _W.end(); ++iPtr)
      std::cout << "\nE = " << iPtr->first * binWidth << ", W = " << iPtr->second;
    std::cout << std::endl;

    //Now navigate to the histogram and load the data
    std::istringstream HistogramData
      (std::string(mainNode.getNode("EnergyHist").getNode("HistogramWeighted")));
    //Load the histogram data from string to array
    histogramEntry tmpData;
    while (HistogramData >> tmpData.X[0])
      {
	for (size_t i = 1; i < NGamma; ++i)
	  HistogramData >> tmpData.X[i];

	HistogramData >> tmpData.Probability;
	
	tmpData.Probability *= binWidth;
	
	data.push_back(tmpData);
      }
  }
  
  bool operator<(const SimulationData& d2) const
  { return gamma[0] < d2.gamma[0]; }

  bool operator>(const SimulationData& d2) const
  { return gamma[0] > d2.gamma[0]; }

  std::string fileName;
  std::vector<long double> gamma;
  long double logZ;
  long double new_logZ;
  long double binWidth;
  bool refZ;

  //Contains the histogram, first axis is bin entry
  //second axis is value of X with the final entry being the probability
  struct histogramEntry
  {
    typedef std::array<long double, NGamma> Xtype;
    Xtype X;
    long double Probability;
  };
  std::vector<histogramEntry> data;

  std::unordered_map<int, double> _W;


  //Bottom and top contain the window of systems used for the calculation
  long double calc_logZ(size_t bottom=0, size_t top=0) const
  {
    long double sum1 = 0.0;
    
    //If top = 0, then use all systems
    if (top == 0) top = SimulationDataData.size() - 1;

    for (size_t i(bottom); i <= top; ++i)
      for (const histogramEntry& simdat : SimulationDataData[i].data)
	{
	  long double sum2 = 0.0;

	  for (size_t j(bottom); j <= top; ++j)
	    {
	      //Determine (\gamma_i- \gamma) \cdot X 
	      long double dot = 0.0;

	      if (NGamma != 1) 
		M_throw() << "For multiple gamma reweighting, one must be designated as E and used in the W lookup";

	      { const int i = 0;
		dot += (SimulationDataData[j].gamma[i] - gamma[i]) * simdat.X[i]
		  + SimulationDataData[j].W(simdat.X[i]) - W(simdat.X[i]);
	      }

	      //This is Z^{-1} \exp[(\gamma_i- \gamma) \cdot X + W_i(X) - W(X)]
	      sum2 += exp(dot - SimulationDataData[j].logZ);
	    }
	  
	  //simdat.Probability is H(X,\gamma), in the input H is
	  //normalised. Thus this assumes that all simulations are of
	  //the same statistical weight! (This is true for results
	  //from a single replica exchange simulation)
	  sum1 += simdat.Probability / sum2;
	}
   
    return log(sum1);
  }
  
  void recalc_newlogZ(size_t bottom, size_t top) 
  { 
    if (!refZ) new_logZ = calc_logZ(bottom, top); 
  }

  long double calc_error()
  { 
    //Return an error of 0 if this is the reference simulation!
    if (refZ) return 0;

    //In case of new_logZ going to zero don't use relative values
    if (new_logZ == 0)
      {
	if (logZ == 0)
	  return 0;
	else
	  return fabs((new_logZ - logZ) / logZ);
      }
    else
      return fabs((new_logZ - logZ) / new_logZ);
  }

  void iterate_logZ() { logZ = new_logZ; }

  inline double W(double E) const 
  { 
    std::unordered_map<int, double>::const_iterator 
      iPtr = _W.find(lrint(E / binWidth));
    if (iPtr != _W.end())
      return -iPtr->second;
    return 0;
  }
};

struct ldbl 
{
  long double val;
  ldbl():val(0) {}
  ldbl(long double nval):val(nval) {}

  operator long double&() { return val; }
  operator const long double&() const { return val; }
};

typedef std::pair<SimulationData::histogramEntry::Xtype, ldbl> densOStatesPair;
typedef std::map<SimulationData::histogramEntry::Xtype, ldbl> densOStatesMap;
typedef std::vector<densOStatesPair> densOStatesType;
densOStatesType densOStates;
  

void
solveWeightsInRange(size_t bottom = 0, size_t top = 0)
{
  //If top = 0, then use all systems
  if (top == 0) top = SimulationDataData.size() - 1;

  double err = 0.0;

  do
    {
      for (size_t i = NStepsPerStep; i != 0; --i)
	{
	  for (size_t i(bottom); i <= top; ++i)
	    SimulationDataData[i].recalc_newlogZ(bottom, top);
	  
	  for (size_t i(bottom); i <= top; ++i)
	    SimulationDataData[i].iterate_logZ();
	}

      //Now the error checking run
      err = 0.0;
      for (size_t i(bottom); i <= top; ++i)
	{
	  SimulationDataData[i].recalc_newlogZ(bottom, top);
	  
	  if (SimulationDataData[i].calc_error() > err)
	    err = SimulationDataData[i].calc_error();
	}

      //May as well use this as an iteration too
      for (size_t i(bottom); i <= top; ++i)
	SimulationDataData[i].iterate_logZ();

      printf("\r%E", err);
      fflush(stdout);
    }
  while(err > minErr);
}

void
solveWeightsPiecemeal()
{
  std::cout << "##################################################\n";
  std::cout << "Solving for Z's, in a rolling piecemeal fashion\n";

  size_t startingPieceSize = 5;

  size_t stoppingPieceSize = SimulationDataData.size() / 2 + 1;
  
  if (startingPieceSize > SimulationDataData.size()) startingPieceSize = SimulationDataData.size();

  for (size_t pieceSize = startingPieceSize; pieceSize < stoppingPieceSize; pieceSize += 5)
    {
      for (SimulationData& simdat : SimulationDataData)
	simdat.refZ = false;

      //Set the first system as the reference point
      SimulationDataData.front().refZ = true;
      
      std::cout << "\rSolving 0 to " << pieceSize << ", Long iteration step\n";
      //Now solve the first piece.
      solveWeightsInRange(0, pieceSize);
      
      for (size_t bottom(1), top(pieceSize + 1); top < SimulationDataData.size(); ++bottom, ++top)
	{
	  //Freeze the lower half of the last set of systems 
	  
	  //This is required for the first iterations but only one system
	  //is frozen per further iteration
	  for (size_t i = bottom-1; i < bottom + (top-bottom)/2; ++i )
	    SimulationDataData[i].refZ = true;
	  
	  std::cout << "\rSolving " << bottom << " to " << top << "\n";
	  solveWeightsInRange(bottom, top);
	}
    }
  
  for (SimulationData& simdat : SimulationDataData)
    simdat.refZ = false;

  //Set the first system as the reference point
  SimulationDataData.front().refZ = true;

  std::cout << "\rFinal Solution step 0 to " << SimulationDataData.size()-1 << "\n";
  solveWeightsInRange();

  std::cout << "\nIteration complete\n";
}


void calcDensityOfStates()
{
  densOStates.clear();

  densOStatesMap accumilator;
  std::cout << "##################################################\n";
  std::cout << "Density of states\n";

  //Sum up the histogram entries for each value of X
  for (const SimulationData& dat : SimulationDataData)
    for (const SimulationData::histogramEntry& simdat : dat.data)
    accumilator[simdat.X] += simdat.Probability;

  //For each X, calculate the density of states
  for (const auto& dat : accumilator)
    {
      //First work out the divisor
      long double sum = 0.0;
      for (const SimulationData& dat2 : SimulationDataData)
	{
	  long double tmp = 0;

	  //Adding in the W factor
	  if (NGamma != 1) 
	    M_throw() << "For multiple gamma reweighting, one must be designated as E and used in the W lookup";

	  tmp += dat2.W(dat.first[0]);

	  for (size_t i(0); i < NGamma; ++i)
	    tmp += dat2.gamma[i] * dat.first[i];

	  sum += exp(tmp - dat2.logZ);
	}
      
      //The sum of the histogram entries is already calculated so just push it on
      densOStates.push_back(std::make_pair(dat.first, dat.second / sum));
    }
}

void outputDensityOfStates()
{
  
  std::fstream of("StateDensity.out",  ios::trunc | ios::out);

  of << std::setprecision(std::numeric_limits<long double>::digits10);

  for (const densOStatesPair& dat : densOStates)
    {
      for (size_t i(0); i < NGamma; ++i)
	of << dat.first[i] << " " ;

      of << dat.second << "\n";
    }
  
  of.close();
}

void outputLogZ()
{
  std::fstream of("logZ.out",  ios::trunc | ios::out);

  of << std::setprecision(std::numeric_limits<long double>::digits10);

  for (const SimulationData& dat : SimulationDataData)
    of << dat.gamma[0] << " " << dat.logZ << "\n";

  of.close();
}

void outputMoments()
{
  std::cout << "##################################################\n";
  std::cout << "Calculating  moments\n";
  
  for (const SimulationData& dat : SimulationDataData)
    {
      std::cout << "Writing " + (dat.fileName + std::string(".ReweightedEnergyHist"))
		<< "\n";

      std::fstream Eof((dat.fileName + std::string(".ReweightedEnergyHist")).c_str(),  
		       ios::trunc | ios::out);

      Eof << std::setprecision(std::numeric_limits<long double>::digits10);
      //Calc Z
      long double Z = 0.0;

      for (const densOStatesPair& dat2 : densOStates)
	{
	  long double tmp(0);
	  for (size_t i(0); i < NGamma; ++i)
	    tmp += dat.gamma[i] * dat2.first[i];

	  Z += std::exp(std::log(dat2.second) + tmp);
	}
      
      Z = log(Z);
      
      //Now calc the normalisation, or first moment
      long double Norm = 0.0;
      for (const densOStatesPair& dat2 : densOStates)
	{
	  long double tmp(0);
	  for (size_t i(0); i < NGamma; ++i)
	    tmp += dat.gamma[i] * dat2.first[i];

	  Norm += std::exp(std::log(dat2.second) + tmp - Z);
	}

      for (const densOStatesPair& dat2 : densOStates)
	{
	  long double tmp(0);

	  for (size_t i(0); i < NGamma; ++i)
	    tmp += dat.gamma[i] * dat2.first[i];

	  for (size_t i(0); i < NGamma; ++i)
	    Eof << dat2.first[i] << " ";

	  Eof << (std::exp(std::log(dat2.second) + tmp - Z) / Norm) / SimulationDataData.front().binWidth << "\n";
	}
    }

  size_t steps = 100 * (SimulationDataData.size()-1) + 1;

  long double stepsize = (betaMax - betaMin) / steps;

  std::vector<std::pair<double,double> > Cv;
  {
    std::fstream Eof("Energy.out",  ios::trunc | ios::out);
    std::fstream E2of("Energy2.out",  ios::trunc | ios::out);
    std::fstream Cvof("Cv.out",  ios::trunc | ios::out);
    
    Eof << std::setprecision(std::numeric_limits<long double>::digits10);
    E2of << std::setprecision(std::numeric_limits<long double>::digits10);
    Cvof << std::setprecision(std::numeric_limits<long double>::digits10);
    
    for (size_t step = 0; step <= steps; step++)
      {
	long double Beta = betaMin + step * stepsize;

	//Calc Z
	long double Z = 0.0;
	for (const densOStatesPair& dat : densOStates)
	  {
	    //ERROR! Not generalised for multipe Histogram arguments!
	    Z += exp (log(dat.second) + Beta * dat.first[0]);
	  }
	
	Z = log(Z);
	
	//Now calc the normalisation, or first moment
	long double Norm = 0.0;
	long double Eavg = 0.0;
	long double E2avg = 0.0;
	
	for (const densOStatesPair& dat : densOStates)
	  {
	    //ERROR! Not generalised for multipe Histogram arguments!
	    long double temp = exp(log(dat.second) + Beta * dat.first[0] - Z);
	    Norm += temp;
	    Eavg += temp * dat.first[0];
	    E2avg += temp * dat.first[0] * dat.first[0];
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
	if ((olddt < 0) != (thisdt < 0))
	  {
	    if (olddt < 0)
	      CvMin << (iPtr-1)->first << " " << (iPtr-1)->second << "\n";
	    else
	      CvMax << (iPtr-1)->first << " " << (iPtr-1)->second << "\n";
	  }
	olddt = thisdt;
      }
    
    CvMax.close();
    CvMin.close();
  }
  
}

int
main(int argc, char *argv[])
{
  std::cout << "dynahist_rw  Copyright (C) 2011  Marcus N Campbell Bannerman\n"
	    << "This program comes with ABSOLUTELY NO WARRANTY.\n"
	    << "This is free software, and you are welcome to redistribute it\n"
	    << "under certain conditions. See the licence you obtained with\n"
	    << "the code\n";

#if !defined(__APPLE__) && !defined(_WIN32)
  //This is so the program crashes out when floating point errors occur
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif

  try {
    namespace po = boost::program_options;
    
    boost::program_options::options_description systemopts("Program Options");
    
    systemopts.add_options()
      ("help", "Produces this message")   
      ("data-file", po::value<std::vector<std::string> >(), "Specify a config file to load, or just list them on the command line")
      ("alpha", po::value<long double>()->default_value(1), "A fraction of the difference between the old and new logZ's to use, use to stop divergence")
      ("NSteps,N", po::value<size_t>()->default_value(10), "Number of steps to take before testing the error and spitting out the current vals")
      ("Tmin", po::value<double>(), "Set the coldest temperature to output calculated data for (Cv.out, Energy.out) etc. If unset this defaults to the temperature of the coldest simulation.")
      ("Tmax", po::value<double>(), "Set the hottest temperature to output calculated data for (Cv.out, Energy.out) etc. If unset this defaults to the temperature of the hottest simulation.")
      ;

    boost::program_options::positional_options_description p;
    p.add("data-file", -1);
    
    boost::program_options::store(po::command_line_parser(argc, argv).
				  options(systemopts).positional(p).run(), vm);
    boost::program_options::notify(vm);
    
    if (vm.count("help") || !vm.count("data-file")) 
      {
	M_throw() << "Usage : dynahist_rw <OPTION>...<data-file(s)>\n"
		  << "Determines the weighting functions for the histograms\n"
		  << systemopts << "\n";
      }

    alpha = vm["alpha"].as<long double>();
    NStepsPerStep = vm["NSteps"].as<size_t>();

    //Data load
    for (std::string fileName : vm["data-file"].as<std::vector<std::string> >())
      SimulationDataData.push_back(SimulationData(fileName));

    for (const SimulationData& dat : SimulationDataData)
      if (dat.binWidth != SimulationDataData.front().binWidth)
	{
	  std::cout << "Not all of the output files have the same bin width for the internal energy histograms!\n Aborting\n\n";
	  return 1;
	}
    //Sort by temperature
    std::sort(SimulationDataData.begin(), SimulationDataData.end());

    //Save the temperature range
    betaMax = SimulationDataData.front().gamma[0];
    if (vm.count("Tmin"))
      {
	if (vm["Tmin"].as<double>() <= 0)
	  M_throw() << "Tmin must be positive and non-zero";
	betaMax = - 1.0 / vm["Tmin"].as<double>();
      }

    betaMin = SimulationDataData.back().gamma[0];
    if (vm.count("Tmin"))
      {
	if (vm["Tmax"].as<double>() <= 0)
	  M_throw() << "Tmin must be positive and non-zero";
	betaMin = - 1.0 / vm["Tmax"].as<double>();
      }

    //Output ordered sims
    std::cout << "##################################################\n";
    for (const SimulationData& dat : SimulationDataData)
      std::cout << dat.fileName << " NData = " << dat.data.size() << " gamma[0] = " << dat.gamma[0] << "\n";

    solveWeightsPiecemeal();
    
    std::cout << "##################################################\n";
    for (const SimulationData& dat : SimulationDataData)
      std::cout << dat.fileName << " logZ = " << dat.logZ << "\n";

    //Now start the fun
    outputLogZ();
    calcDensityOfStates();
    outputDensityOfStates();
    outputMoments();
    
  }
  catch (std::exception& cep)
    {
      fflush(stdout);
      std::cerr << cep.what() << "\nMAIN: Reached Main Error Loop\n";
      return 1;
    }
}
