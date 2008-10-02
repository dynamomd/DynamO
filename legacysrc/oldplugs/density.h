#ifndef COPDENS_H
#define COPDENS_H

#include "outputplugin.h"
#include "../extcode/xmlwriter.h"

#define bincount 200

class COPDensity: public COutputPlugin
{
 public:
  COPDensity(const std::vector<CParticle> &, const CDynamics * const);
  ~COPDensity();

  void collisionUpdate(const CCollision &, const CEventData &b);

  void output(xmlw::XmlStream &); 

  void periodicOutput();

  virtual COutputPlugin *Clone() { return new COPDensity(*this); };
  
 protected:

  CVector<long> bin[bincount];
  double binwidth;
  long samplecount;
};

#endif
