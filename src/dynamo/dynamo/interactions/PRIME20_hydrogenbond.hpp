#pragma once

#include <dynamo/interactions/captures.hpp>
#include <dynamo/interactions/glyphrepresentation.hpp>
#include <dynamo/simulation.hpp>

namespace dynamo {
   class IPRIME20_HydrogenBond: public ISingleCapture, public GlyphRepresentation
   {
   public:
      template<class T>
      IPRIME20_HydrogenBond(dynamo::Simulation* tmp, T wd, C2Range* nR, std::string name):
         ISingleCapture(tmp,nR),
         _wellDepth (Sim->_properties.getProperty (wd, Property::Units::Energy        ( ))),
      { intName = name; }

      IPRIME20_HydrogenBond(const magnet::xml::Node&, dynamo::Simulation*);

      void operator<<(const magnet::xml::Node&);

      virtual size_t glyphsPerParticle() const { return 1; }
      virtual Vector getGlyphSize(size_t ID, size_t subID) const;
      virtual Vector getGlyphPosition(size_t ID, size_t subID) const;

      virtual double getExcludedVolume(size_t) const;

      virtual double maxIntDist() const;

      virtual void checkOverlaps(const Particle&, const Particle&) const;

      virtual bool captureTest(const Particle&, const Particle&) const;

      virtual void initialise(size_t);

      virtual IntEvent getEvent(const Particle&, const Particle&) const;

      virtual void runEvent(Particle&, Particle&, const IntEvent&) const;

      virtual void outputXML(magnet::xml::XmlStream&) const;

      virtual double getInternalEnergy() const;

      virtual double getInternalEnergy(const Particle&, const Particle&) const;

   protected:
      IPRIME20_HydrogenBond(dynamo::Simulation* tmp, C2Range* nR):
         ISingleCapture(tmp,nR) {}

      shared_ptr<Property> _wellDepth;
      shared_ptr<Property> _diameter; //TODO not sure if this supports different values for different pairs in the the 6 atom set

      //Function to return array of relevant IDs:
      virtual void ID_array(const Particle&, const Particle&, Particle&[]);
   };
}
