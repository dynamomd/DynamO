#include "RenderObj.hpp"

RenderObj::RenderObj(bool hostTransfers):
  _RenderMode(TRIANGLES),
  _hostTransfers(hostTransfers)
{}

RenderObj::~RenderObj()
{}

