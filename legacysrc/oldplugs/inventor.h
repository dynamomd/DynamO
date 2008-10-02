#ifndef COPINVENTOR_H
#define COPINVENTOR_H

#include "outputplugin.h"
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/Xt/SoXt.h>

class COPInventor: public COutputPlugin
{
 public:
  COPInventor(const std::vector<CParticle> &, const CDynamics * const);
  ~COPInventor();
  
  void collisionUpdate(const CCollision &, const CEventData &);

  virtual COutputPlugin *Clone() { return new COPInventor(*this); };

 private:
  int frameCount;

  SoSeparator *makeMoleculeScene();
  void viewScene();
  Widget myWindow;
};

#endif
