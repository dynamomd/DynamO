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

#ifndef CEReplexer_H
#define CEReplexer_H

#include "engine.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>

class CEReplexer: public CEngine
{
public:
  CEReplexer(const boost::program_options::variables_map&);

  virtual ~CEReplexer() {}

  virtual void printStatus();

  virtual void runSimulation();

  virtual void initialisation();
  
  virtual void forceShutdown();

  virtual void peekData();

  virtual void finaliseRun() {}

  static void getOptions(boost::program_options::options_description&);

  virtual void outputData();
  
  virtual void outputConfigs();

protected:

  typedef enum {
    NoSwapping = 0,
    AlternatingSequence = 1,
    SinglePair = 2,
    RandomPairs = 3,
    RandomSelection = 4
  } Replex_Mode_Type;

  struct simData
  {
    explicit simData(int ID, Iflt rT):
      simID(ID), swaps(0), attempts(0), upSims(0), downSims(0),
      realTemperature(rT)
    {}
    bool operator<(const simData& sdat) const
    {return simID < sdat.simID; }
    
    int simID;
    size_t swaps;
    size_t attempts;
    size_t upSims;
    size_t downSims;
    Iflt realTemperature;
  };

  typedef std::pair<Iflt, simData> replexPair;

  boost::scoped_array<CSimulation> Simulations;
  Iflt replicaEndTime;  
  Replex_Mode_Type ReplexMode;
  std::vector<replexPair> temperatureList;
  std::vector<int> SimDirection;
  std::vector<char> roundtrip;
  size_t replexSwapCalls;
  size_t round_trips;
  boost::posix_time::ptime start_Time;
  boost::posix_time::ptime end_Time;
  bool SeqSelect;
  unsigned int nSims;
  bool peekMode;

  virtual void preSimInit();

  virtual void setupSim(CSimulation&, const std::string);

  virtual void postSimInit(CSimulation&);

  //Replica Exchange code
  void ReplexInit();
  void ReplexSwap(Replex_Mode_Type localMode);
  void ReplexConfigOutput(std::vector<std::string>&);
  void ReplexDataOutput(std::vector<std::string>&);
  void ReplexSwapTicker();
  void AttemptSwap(unsigned int, unsigned int);
  inline void AttemptSwap(unsigned int ID) { AttemptSwap(ID, ID + 1); }
};

#endif
