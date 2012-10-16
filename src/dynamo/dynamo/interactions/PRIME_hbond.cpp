//#include <dynamo/interactions/PRIME_hbond.hpp>
//#include <dynamo/BC/BC.hpp>
//
//#include <dynamo/units/units.hpp>
//#include <dynamo/globals/global.hpp>
//#include <dynamo/particle.hpp>
//#include <dynamo/interactions/intEvent.hpp>
//#include <dynamo/species/species.hpp>
//#include <dynamo/2particleEventData.hpp>
//#include <dynamo/dynamics/dynamics.hpp>
//#include <dynamo/simulation.hpp>
//#include <dynamo/schedulers/scheduler.hpp>
//#include <dynamo/NparticleEventData.hpp>
//#include <dynamo/outputplugins/outputplugin.hpp>
//#include <magnet/xmlwriter.hpp>
//#include <magnet/xmlreader.hpp>
//#include <cmath>
//#include <iomanip>
//
//namespace dynamo {
//
//   //////////////////////////////////////////////////
//   //                Initialisation                //
//   //////////////////////////////////////////////////
//
//   IPRIME_Hbond::IPRIME_Hbond(const magnet::xml::Node& XML, dynamo::Simulation* tmp): ISingleCapture(tmp, NULL) {
//         operator<<(XML);                       // A temporary value!
//         ID_array(temp1, temp2, this->IDs, this->IDpairs); // Grab arrays of IDs.
//   }
//
//   void IPRIME_Hbond::operator<<(const magnet::xml::Node& XML) {
//      if (strcmp(XML.getAttribute("Type"),"PRIME_Hbond"))
//         M_throw() << "Attempting to load PRIME_Hbond from non PRIME_Hbond entry";
//
//      Interaction::operator<<(XML);
//
//      try {
//         _wellDepth = Sim->_properties.getProperty(XML.getAttribute("WellDepth"), Property::Units::Energy());
//
//         intName = XML.getAttribute("Name");
//         ISingleCapture::loadCaptureMap(XML);
//      }
//      catch (boost::bad_lexical_cast &) { M_throw() << "Failed a lexical cast in IPRIME_Hbond"; }
//   }
//
//   void IPRIME_Hbond::initialise(size_t nID) { ID = nID; ISingleCapture::initCaptureMap(); }
//
//   //////////////////////////////////////////////////
//   //               Single-ID methods              //
//   //////////////////////////////////////////////////
//
//   Vector IPRIME_Hbond::getGlyphSize(size_t ID, size_t subID) const {
//      double diam = _diameter->getProperty(ID);
//      return Vector(diam, diam, diam)
//   }
//
//   Vector IPRIME_Hbond::getGlyphPosition(size_t ID, size_t subID) const {
//      Vector retval = Sim->particles[ID].getPosition();
//      Sim->BCs->applyBC(retval);
//      return retval;
//   }
//
//   double IPRIME_Hbond::getExcludedVolume(size_t ID) const {
//      double diam = _diameter->getProperty(ID);
//      return diam * diam * diam * M_PI / 6.0; 
//   }
//
//   ////////////////////////////////////////////////
//   //             Interaction methods            //
//   ////////////////////////////////////////////////
//
//   //Not sure if this is any use:
//   double IPRIME_Hbond::maxIntDist() const { return _diameter->getMaxValue() * _lambda->getMaxValue(); }
//
//   virtual void ID_array(const Particle &temp1, const Particle &temp2, Particle &IDs[6], Particle &IDpairs[10]){
//      //Writes the particles to the array in order: CO, CO's CH, CO's NH, NH, NH's CO, NH's CH.
//      Particle &CO, &NH, &CO_bond[2], &NH_bond[2];
//
//      //Identify CO and NH
//      (<temp1.type> == "CO") ? ( IDs[0] = temp1 && IDs[3] = temp2 ) : ( IDs[3] = temp1 && IDs[0] = temp2 );
//
//      //Identify CO's bonded neighs
//      ( <IDs[0].bondedneighs[0].type> == "CH" ) ? ( IDs[1] = IDs[0].bondedneighs[0] && IDs[2] = IDs[0].bondedneighs[1] )
//         : ( IDs[2] = IDs[0].bondedneighs[0] && IDs[1] = IDs[0].bondedneighs[1] )
//
//      //Identify NH's bonded neighs
//      ( <IDs[3].bondedneighs[0].type> == "CO" ) ? ( IDs[4] = IDs[3].bondedneighs[0] && IDs[5] = IDs[3].bondedneighs[1] )
//         : ( IDs[5] = IDs[3].bondedneighs[0] && IDs[4] = IDs[3].bondedneighs[1] )
//
//      //Interaction pairs
//      IDpairs[0] = IDs[0]; IDpairs[1] = IDs[3]; // CO and NH
//      IDpairs[2] = IDs[0]; IDpairs[3] = IDs[4]; // CO and NH's CO
//      IDpairs[4] = IDs[0]; IDpairs[5] = IDs[5]; // CO and NH's CH
//      IDpairs[6] = IDs[3]; IDpairs[7] = IDs[1]; // NH and CO's CH
//      IDpairs[8] = IDs[3]; IDpairs[9] = IDs[2]; // NH and CO's NH
//
//      //TODO add some debug stuff: check for all types non-lazily, make sure no extra bonds, etc.
//   }
//
//   bool IPRIME_Hbond::captureTest(const Particle& temp1, const Particle& temp2) const {
//      if (&(*(Sim->getInteraction(temp1, temp2))) != this) return false;
//
//      double d, l;
//      bool captured = true;
//
//      //5 constraints to satisfy
//      for (int i = 0; i < 10 && captured; i+=2){
//         d = ( _diameter->getProperty(this->IDpairs[i].getID()) + _diameter->getProperty(this->IDpairs[i+1].getID()) ) * 0.5;
//         l = ( _lambda->getProperty(this->IDpairs[i].getID()) + _lambda->getProperty(this->IDpairs[i+1].getID()) ) * 0.5;
//
//         //Only captured if all  5 constraints are met hence ⋀
//         captured = captured && Sim->dynamics->sphereOverlap(this->IDpairs[i], this->IDpairs[i+1], d * l);
//
//#ifdef DYNAMO_DEBUG
//         if (Sim->dynamics->sphereOverlap(this->IDpairs[i], this->IDpairs[i+1], d)){
//            derr << "Warning! Two particles might be overlapping" << "Overlap is " << Sim->dynamics->sphereOverlap(this->IDpairs[i], this->IDpairs[i+1], d) 
//      / Sim->units.unitLength() << "\nd = " << d / Sim->units.unitLength() << std::endl;
//         }
//#endif
//      }
//      return captured;
//   }
//
//   IntEvent IPRIME_Hbond::getEvent(const Particle &p1, const Particle &p2) const {
//      //Interaction pairs are: CO & NH, CO & NH_bond[0], CO & NH_bond[0], NH & CO_bond[0], NH & CO_bond[1]
//#ifdef DYNAMO_DEBUG
//      if (!Sim->dynamics->isUpToDate(p1)) M_throw() << "particle 1 is not up to date";
//
//      if (!Sim->dynamics->isUpToDate(p2)) M_throw() << "particle 2 is not up to date";
//
//      if (p1 == p2) M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
//#endif 
//
//      double d,l;
//      double dt_IN = 0, dt_OUT = HUGE_VAL, dt_CORE, dt_temp; //Time until next IN event.
//      bool = captured;
//      IntEvent OUT_event, IN_event, CORE_event;
//
//      //Cycle through all pairs and look for important events. Events occur when:
//      //1 - CO and NH have a CORE event. 2 - Any particle has a WELL_OUT event IF all particles were previously WELL_IN
//      //3 - Any particle has a WELL_IN event while all other particles are WELL_IN.
//      //Need to return the soonest of these 3 possibilities.
//
//      for (int i = 0; i < 10; i+=2){
//         d = (_diameter->getProperty(this->IDpairs[i].getID()) + _diameter->getProperty(this->IDpairs[i+1].getID())) * 0.5;
//         l = (_lambda->getProperty(this->IDpairs[i].getID()) + _lambda->getProperty(this->IDpairs[i+1].getID())) * 0.5;
//         captured = isCaptured(this->IDpairs[i], this->IDpairs[i+1]);
//
//         //Look for making or breaking the full h-bond:
//         if (uncaptured == 0 && captured) {
//            dt_temp = Sim->dynamics->SphereSphereOutRoot(this->IDpairs[i], this->IDpairs[i+1], l * d);
//            if (dt_temp < dt_OUT){ // Only care about first OUT event.
//               dt_OUT    = dt_temp;
//               OUT_event = IntEvent(this->IDpairs[i], this->IDpairs[i+1], dt_OUT, WELL_OUT, *this);
//            }
//         } else if (!captured) {
//            uncaptured++;
//            dt_temp = Sim->dynamics->SphereSphereInRoot(this->IDpairs[i], this->IDpairs[i+1], l * d);
//            if (dt_temp > dt_IN){ // Only care about last IN event.
//               dt_IN    = dt_temp;
//               IN_event = IntEvent(this->IDpairs[i], this->IDpairs[i+1], dt, WELL_IN, *this);
//      }  }  }
//
//      dt_CORE    = Sim->dynamics->SphereSphereInRoot(p1, p2, d);
//      CORE_event = IntEvent(this->IDpairs[i], this->IDpairs[i+1], dt_CORE, CORE, *this);
//
//      if (uncaptured == 0){ //IN right now; event could be OUT or CORE
//         return (dt_IN < dt_CORE) ? IN_event : CORE_event;
//      }
//
//      if (uncaptured != 0){ //OUT right now; event could be IN or CORE
//         return (dt_OUT < dt_CORE) ? OUT_event : CORE_event;
//   }  }
//
//   void IPRIME_Hbond::checkOverlaps(const Particle& p1, const Particle& p2) const {
//      //sites other than the main two will be done via their own 'HardSphere' interactions; not here.
//
//      Vector  rij = part1.getPosition() - part2.getPosition();
//      Sim->BCs->applyBC(rij);
//      double r2 = rij.nrm2();
//
//      double d = (_diameter->getProperty(part1.getID()) + _diameter->getProperty(part2.getID())) * 0.5;
//      double l = (_lambda->getProperty(part1.getID()) + _lambda->getProperty(part2.getID())) * 0.5;
//
//      double d2 = d * d;
//      double ld2 = d2 * l * l;
//
//      if (isCaptured(part1, part2)) {
//         if (r2 < d2)
//            derr << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID()
//            << ", ID2=" << part2.getID() << "\nR_ij^2=" << r2 / pow(Sim->units.unitLength(),2) << "\nd^2=" 
//            << d2 / pow(Sim->units.unitLength(),2) << std::endl;
//
//         if (r2 > ld2)
//            derr << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID()
//            << ", ID2=" << part2.getID() << "\nR_ij^2=" << r2 / pow(Sim->units.unitLength(),2)
//            << "\n(lambda * d)^2=" << ld2 / pow(Sim->units.unitLength(),2) << std::endl;
//      }  else if (r2 < ld2) {
//         if (r2 < d2)
//            derr << "Overlap error\n ID1=" << part1.getID() << ", ID2=" << part2.getID() << "\nR_ij^2="
//            << r2 / pow(Sim->units.unitLength(),2) << "\n(d)^2=" << d2 / pow(Sim->units.unitLength(),2) << std::endl;
//         else
//            derr << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
//            << ", ID2=" << part2.getID() << "\nR_ij^2=" << r2 / pow(Sim->units.unitLength(),2)
//            << "\n(lambda * d)^2=" << ld2 / pow(Sim->units.unitLength(),2) << std::endl;
//      }
//   }
//
//   void IPRIME_Hbond::outputXML(magnet::xml::XmlStream& XML) const {
//      XML << magnet::xml::attr("Type") << "PRIME_Hbond"
//      << magnet::xml::attr("WellDepth") << _wellDepth->getName()
//      << magnet::xml::attr("Name") << intName
//      << *range;
//      ISingleCapture::outputCaptureMap(XML);
//   }
//
//   double IPRIME_Hbond::getInternalEnergy(const Particle& p1, const Particle& p2) const {
//      return - 0.5 * (_wellDepth->getProperty(p1.getID()) +_wellDepth->getProperty(p2.getID())) * isCaptured(p1, p2);
//   }
//   ////////////////////////////////
//   //             WIP            //
//   ////////////////////////////////
//
//   void ISquareWell::runEvent(Particle& p1, Particle& p2, const IntEvent& iEvent) const {
//      ++Sim->eventCount;
//
//      switch( iEvent.getType() ){
//
//         case CORE: { //CORE events can only be between the NH and CO.
//               PairEventData retVal(Sim->dynamics->SmoothSpheresColl(iEvent, 1.0, d2, CORE));
//               Sim->signalParticleUpdate(retVal);
//
//               Sim->ptrScheduler->fullUpdate(p1, p2);
//
//               BOOST_FOREACH(shared_ptr<OutputPlugin> & Ptr, Sim->outputPlugins)
//                 Ptr->eventUpdate(iEvent, retVal);
//
//               break;
//            }
//
//         case WELL_IN: { //Could be any of the interacting pairs. Query iEvent to find out I guess.
//               //do things
//               break;
//            }
//
//         case WELL_OUT: {
//               //do things
//               break;
//            }
//
//         default:
//            M_throw() << "Unknown collision type";
//      }
//   }
//
//   double ISquareWell::getInternalEnergy() const {
//      //TODO not sure how to edit this.
//      //Once the capture maps are loaded just iterate through that determining energies
//      double Energy = 0.0;
//      typedef std::pair<size_t, size_t> locpair;
//
//      BOOST_FOREACH(const locpair& IDs, captureMap)
//         Energy += 0.5 * (_wellDepth->getProperty(IDs.first) +_wellDepth->getProperty(IDs.second));
//
//      return -Energy;
//   }
//
//
//}
