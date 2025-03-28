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
#include <dynamo/BC/BC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/trajectory.hpp>
#include <dynamo/systems/system.hpp>
#include <dynamo/units/units.hpp>
#include <iomanip>

namespace dynamo {
OPTrajectory::OPTrajectory(const dynamo::Simulation *t1,
                           const magnet::xml::Node &)
    : OutputPlugin(t1, "Trajectory") {}

OPTrajectory::OPTrajectory(const OPTrajectory &trj) : OutputPlugin(trj) {
  if (trj.logfile.is_open())
    trj.logfile.close();
}

void OPTrajectory::initialise() {
  if (logfile.is_open())
    logfile.close();

  logfile.open("trajectory.out", std::ios::out | std::ios::trunc);
  logfile.precision(4);
  logfile.setf(std::ios::fixed);
}

void OPTrajectory::eventUpdate(const Event &eevent, const NEventData &SDat) {
  logfile << std::setw(8) << std::setfill('0') << Sim->eventCount
          << ", Source=" << eevent._source << ", SourceID=" << eevent._sourceID
          << ", Event Type=" << eevent._type
          << ", t=" << Sim->systemTime / Sim->units.unitTime()
          << ", dt=" << eevent._dt / Sim->units.unitTime();

  for (const ParticleEventData &pData : SDat.L1partChanges) {
    logfile << "\n";
    const Particle &part = Sim->particles[pData.getParticleID()];
    logfile << "   1PEvent: p1=" << part.getID()
            << ", Type=" << pData.getType();
    Vector delP = Sim->species[pData.getSpeciesID()]->getMass(part.getID()) *
                  (part.getVelocity() - pData.getOldVel());
    delP /= Sim->units.unitMomentum();
    Vector pos = part.getPosition() / Sim->units.unitLength();
    Vector oldv = pData.getOldVel() / Sim->units.unitVelocity();
    Vector newv = part.getVelocity() / Sim->units.unitVelocity();
    logfile << ", delP1=" << delP.toString() << ", pos=" << pos.toString()
            << ", vel=" << newv.toString() << ", oldvel=" << oldv.toString()
            << "\n";
  }

  for (const PairEventData &pData : SDat.L2partChanges) {
    const size_t id1 = std::min(pData.particle1_.getParticleID(),
                                pData.particle2_.getParticleID());
    const size_t id2 = std::max(pData.particle1_.getParticleID(),
                                pData.particle2_.getParticleID());
    Vector rij = Sim->particles[id1].getPosition() -
                 Sim->particles[id2].getPosition(),
           vij = Sim->particles[id1].getVelocity() -
                 Sim->particles[id2].getVelocity();

    Sim->BCs->applyBC(rij, vij);
    rij /= Sim->units.unitLength();
    vij /= Sim->units.unitVelocity();

    logfile << "\n   2PEvent:";
    logfile << " p1=" << std::setw(5) << id1 << ", p2=" << std::setw(5) << id2
            << ", delP1="
            << ((id1 == pData.particle1_.getParticleID())
                    ? pData.impulse.toString()
                    : (-pData.impulse).toString())
            << ", |r12|=" << std::setw(5) << rij.nrm()
            << ", post-r12=" << rij.toString()
            << ", post-v12=" << vij.toString()
            << ", post-rvdot=" << (vij | rij);
  }
  logfile << "\n";
}

void OPTrajectory::output(magnet::xml::XmlStream &XML) {}
} // namespace dynamo
