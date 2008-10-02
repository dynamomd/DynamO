#ifndef COPCLUSTER_H
#define COPCLUSTER_H

#include "outputplugin.h"
#include "../extcode/xmlwriter.h"

#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/Xt/SoXt.h>

class COPCluster: public COutputPlugin
{
 public:
  COPCluster(const std::vector<CParticle> &, const CDynamics * const );
  ~COPCluster();
 
  void collisionUpdate(const CCollision &, const CEventData &b);

  void output(xmlw::XmlStream &); 

  void periodicOutput();

  SoSeparator *makeScene(double length);

  void viewClusters(double length);

  void writeToFile(SoSeparator *, char *);

  double orderParameter(double);

  virtual COutputPlugin *Clone() { return new COPCluster(*this); };
  
 protected:
};

#endif
