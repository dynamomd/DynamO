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

#include <dynamo/outputplugins/tickerproperty/craig.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <magnet/xmlwriter.hpp>
#include <fstream>
#include <sstream>

// Creates a temperature profile of particles along the x-axis of the
// simulations. This plugin calculates a rolling of the temperature in a
// specified number of bins averaged at every time step.

namespace dynamo {
      OPCraig::OPCraig(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
            OPTicker(tmp,"Craig"),
            nBins(100),
            tickCount(0)
      {
            operator<<(XML);
      }

      void
      OPCraig::operator<<(const magnet::xml::Node& XML)
      {
            //Here is where you process options for the output plugin
            try {
                  if (XML.hasAttribute("numberOfBins")) {
                        nBins = XML.getAttribute("numberOfBins").as<double>();
                  }
            }
            catch (std::exception& e) {
                  M_throw() << "Error while parsing output plugin options\n" << e.what();
            }
      }

      void
      OPCraig::initialise()
      {
            temperatures.resize(nBins);
            densities.resize(nBins);
            for (size_t i = 0; i < nBins; i++) {
                  temperatures[i] = 0.0;
                  densities[i] = 0.0;
            }
            //This is called once, after the simulation is set up, just before
            //the first event is run.
            ticker();
      }

      void
      OPCraig::ticker()
      {
            tickCount++;
            std::vector<double> currentTemperature(nBins, 0.0);
            std::vector<double> currentDensity(nBins, 0.0);
            for (const Particle& p : Sim->particles) {
                  size_t binNumber = floor((0.5 + p.getPosition()[X] / Sim->primaryCellSize[X]) * nBins);
                  if (binNumber > 100) {
                        binnumbers.push_back(p.getPosition()[X]);
                        binnumbers.push_back(Sim->primaryCellSize[X] / 2);
                  }
                  //currentTemperature[binNumber] = 1;
                  // try {
                  //       currentTemperature[binNumber] = 1;
                  // }
                  // catch (std::exception& e) {
                  //       M_throw() << "Error while adding temperature\n" << e.what();
                  // }
            }
            for (size_t i = 0; i < nBins; i++) {
                  //temperatures[i] += currentTemperature[i];
            }
            //This is called periodically, as set by the -t option of dynarun
      }

      void
      OPCraig::output(magnet::xml::XmlStream& XML)
      {
            XML << magnet::xml::tag("Profiles")
                << magnet::xml::attr("NumberOfBins")
                << nBins
                << magnet::xml::attr("BinWidth")
                << Sim->primaryCellSize[0] / nBins
                << magnet::xml::attr("floor_test")
                << floor(2.9);

            XML << magnet::xml::tag("binNumbers")
                << magnet::xml::chardata();
            for (size_t i = 0; i < binnumbers.size(); i++) {
                  XML << binnumbers[i] << " ";
            }
            XML << magnet::xml::endtag("binNumbers");

            XML << magnet::xml::tag("Temperature")
                << magnet::xml::chardata();
            for (size_t i = 0; i < nBins; i++) {
                  XML << temperatures[i] << " ";
            }
            XML << magnet::xml::endtag("Temperature");

            XML << magnet::xml::tag("Density")
                << magnet::xml::chardata();
            for (size_t i = 0; i < nBins; i++) {
                  XML << densities[i] << " ";
            }
            XML << magnet::xml::endtag("Density");

            XML << magnet::xml::endtag("Profiles");
      }

      double
      OPCraig::getTemperature(const Vector& velocity)
      {
            return 0;
      }
}
