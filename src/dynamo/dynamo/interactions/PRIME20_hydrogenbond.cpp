//TODO: make a methods that returns all 6 IDs in a consistent order.
#include <dynamo/interactions/PRIME20_hydrogenbonds.hpp>
#include <dynamo/BC/BC.hpp>

#include <dynamo/units/units.hpp>
#include <dynamo/globals/global.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>
#include <iomanip>

namespace dynamo {

   //################################################
   //#                Initialisation                #
   //################################################

   IPRIME20_HydrogenBond::IPRIME20_HydrogenBond(const magnet::xml::Node& XML, dynamo::Simulation* tmp):
      ISingleCapture(tmp, NULL) { operator<<(XML); } //A temporary value!

   void IPRIME20_HydrogenBond::operator<<(const magnet::xml::Node& XML) {
      if (strcmp(XML.getAttribute("Type"),"PRIME20_HydrogenBond"))
         M_throw() << "Attempting to load PRIME20_HydrogenBond from non PRIME20_HydrogenBond entry";

      Interaction::operator<<(XML);

      try {
         _wellDepth = Sim->_properties.getProperty(XML.getAttribute("WellDepth"), Property::Units::Energy());

         intName = XML.getAttribute("Name");
         ISingleCapture::loadCaptureMap(XML);
      }
      catch (boost::bad_lexical_cast &)
         { M_throw() << "Failed a lexical cast in IPRIME20_HydrogenBond"; }
   }

   void IPRIME20_HydrogenBond::initialise(size_t nID) {
      ID = nID;
      ISingleCapture::initCaptureMap();
  }

   //################################################
   //#               Single-ID methods              #
   //################################################

   Vector IPRIME20_HydrogenBond::getGlyphSize(size_t ID, size_t subID) const {
      double diam = _diameter->getProperty(ID);
      return Vector(diam, diam, diam)
   }

   Vector IPRIME20_HydrogenBond::getGlyphPosition(size_t ID, size_t subID) const {
      Vector retval = Sim->particles[ID].getPosition();
      Sim->BCs->applyBC(retval);
      return retval;
   }

   double IPRIME20_HydrogenBond::getExcludedVolume(size_t ID) const {
      double diam = _diameter->getProperty(ID);
      return diam * diam * diam * M_PI / 6.0; 
   }

   //##############################################
   //#             Interaction methods            #
   //##############################################

   double IPRIME20_HydrogenBond::maxIntDist() const { return _diameter->getMaxValue() * _lambda->getMaxValue(); }

   bool IPRIME20_HydrogenBond::captureTest(const Particle& temp1, const Particle& temp2) const {
      if (&(*(Sim->getInteraction(temp1, temp2))) != this) return false;

      const Particle &CO, &NH, &CO_bond[2], &NH_bond[2];

      //Identify CO NH, and their bonded two neighbors
      (<name of species type of temp1> == "CO") ? ( CO = temp1 && NH = temp2 ) : ( NH = temp1 && CO = temp2 );
      NH_bond[0] = <bonded neighbours of NH>[0];
      NH_bond[1] = <bonded neighbours of NH>[1];
      CO_bond[0] = <bonded neighbours of CO>[0];
      CO_bond[1] = <bonded neighbours of CO>[1];

      const Particle &interaction_pairs[] = {CO, NH, CO, NH_bond[0], CO, NH_bond[0], NH, CO_bond[0], NH, CO_bond[1]};

      //5 constraints to satisfy
      double d, l, captureTest;
      for (int i = 0; i < 10; i+=2){
         d = ( _diameter->getProperty(pairs[i].getID()) + _diameter->getProperty(pairs[i+1].getID()) ) * 0.5;
         l = ( _lambda->getProperty(pairs[i].getID()) + _lambda->getProperty(pairs[i+1].getID()) ) * 0.5;

         captureTest *= Sim->dynamics->sphereOverlap(pairs[i], pairs[i+1], d * l);

#ifdef DYNAMO_DEBUG
         if (Sim->dynamics->sphereOverlap(p1, p2, d)){
            derr << "Warning! Two particles might be overlapping" << "Overlap is " << Sim->dynamics->sphereOverlap(p1, p2, d) 
      / Sim->units.unitLength() << "\nd = " << d / Sim->units.unitLength() << std::endl;
         }
#endif

      }

      return captureTest;
   }

   //##############################
   //#             WIP            #
   //##############################

   IntEvent ISquareWell::getEvent(const Particle &p1, const Particle &p2) const 
   {
#ifdef DYNAMO_DEBUG
      if (!Sim->dynamics->isUpToDate(p1))
         M_throw() << "Particle 1 is not up to date";

      if (!Sim->dynamics->isUpToDate(p2))
         M_throw() << "Particle 2 is not up to date";

      if (p1 == p2)
         M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

      double d = (_diameter->getProperty(p1.getID())
      + _diameter->getProperty(p2.getID())) * 0.5;

      double l = (_lambda->getProperty(p1.getID())
      + _lambda->getProperty(p2.getID())) * 0.5;

      IntEvent retval(p1, p2, HUGE_VAL, NONE, *this);

      if (isCaptured(p1, p2))
         {
   double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
   if (dt != HUGE_VAL)
     {
#ifdef DYNAMO_OverlapTesting
   if (Sim->dynamics->sphereOverlap(p1, p2, d))
     M_throw() << "Overlapping particles found"
          << ", particle1 " << p1.getID()
          << ", particle2 " << p2.getID()
          << "\nOverlap = " 
          << Sim->dynamics.getDynamics()
       .sphereOverlap(p1, p2, d)
       / Sim->units.unitLength();
#endif       
       retval = IntEvent(p1, p2, dt, CORE, *this);
     }

   dt = Sim->dynamics->SphereSphereOutRoot(p1, p2, l * d);
   if (retval.getdt() > dt)
       retval = IntEvent(p1, p2, dt, WELL_OUT, *this);
         }
      else
         {
   double dt = Sim->dynamics->SphereSphereInRoot(p1, p2, l * d);

         if (dt != HUGE_VAL)
   {
#ifdef DYNAMO_OverlapTesting
     if (Sim->dynamics->sphereOverlap(p1, p2, l * d))
       {
         if (Sim->dynamics->sphereOverlap(p1, p2, d))
      M_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
           << "\nparticle1 " << p1.getID() << ", particle2 " 
           << p2.getID() << "\nOverlap = " 
           << Sim->dynamics->sphereOverlap(p1, p2, d)
        / Sim->units.unitLength();
         else
      M_throw() << "Overlapping wells (but not registerd as captured) particles found" 
           << "\nparticle1 " << p1.getID() << ", particle2 " 
           << p2.getID() << "\nOverlap = " 
           << Sim->dynamics->sphereOverlap(p1, p2, l * d)
        / Sim->units.unitLength();

       }
#endif
     retval = IntEvent(p1, p2, dt, WELL_IN, *this);
   }
         }

      return retval;
   }
}
