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
            //This is called once, after the simulation is set up, just before
            //the first event is run.
            temperatures.resize(nBins);
            densities.resize(nBins);
            for (size_t i = 0; i < nBins; i++) {
                  temperatures[i] = 0.0;
                  densities[i] = 0.0;
            }
            ticker();
      }

      void
      OPCraig::ticker()
      {
            //This is called periodically, as set by the -t option of dynarun
            tickCount++;
            for (const Particle& p : Sim->particles) {
                  Vector pos = p.getPosition();
                  Sim->BCs->applyBC(pos);
                  size_t binNumber = floor((0.5 + pos[X] / Sim->primaryCellSize[X]) * nBins);
                  densities[binNumber] += 1.0;
                  temperatures[binNumber] += getTemperature(p.getVelocity(), Sim->species(p)->getMass(p.getID()));
            }
      }

      void
      OPCraig::output(magnet::xml::XmlStream& XML)
      {
            std::vector<double> outputTemperature = temperatures;
            std::vector<double> outputDensity = densities;

            for (size_t i = 0; i < nBins; i++) {
                  if (outputTemperature[i] > 0) {
                        outputTemperature[i] /= 3.0 * outputDensity[i];
                  }
                  outputDensity[i] /= tickCount;
                  outputDensity[i] *= nBins / (volume(Sim->primaryCellSize));
            }

            // GUF
            for (const Particle& p : Sim->particles) {
                  guf.push_back(Sim->species(p)->isSpecies(p));
                  //temperatures[binNumber] += getTemperature(p.getVelocity(), Sim->species(p)->getMass(p.getID()));
            }

            XML << magnet::xml::tag("Profiles")
                << magnet::xml::attr("NumberOfBins")
                << nBins
                << magnet::xml::attr("BinWidth")
                << Sim->primaryCellSize[0] / nBins;

            XML << magnet::xml::tag("Guf")
                << magnet::xml::chardata();
            for (size_t i = 0; i < guf.size(); i++) {
                  XML << guf[i] << " ";
            }
            XML << magnet::xml::endtag("Guf");


            XML << magnet::xml::tag("Temperature")
                << magnet::xml::chardata();
            for (size_t i = 0; i < nBins; i++) {
                  XML << outputTemperature[i] << " ";
            }
            XML << magnet::xml::endtag("Temperature");

            XML << magnet::xml::tag("Density")
                << magnet::xml::chardata();
            for (size_t i = 0; i < nBins; i++) {
                  XML << outputDensity[i] << " ";
            }
            XML << magnet::xml::endtag("Density");

            XML << magnet::xml::endtag("Profiles");
      }

      double
      OPCraig::getTemperature(const Vector& velocity, const double mass)
      {
            return mass * (pow(velocity[0], 2) + pow(velocity[1], 2) + pow(velocity[2], 2));
      }

      double
      OPCraig::volume(const Vector& simulationLength)
      {
            return (simulationLength[0] * simulationLength[1] * simulationLength[2]);
      }
}
