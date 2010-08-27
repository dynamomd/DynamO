#pragma once
#include "RenderObj.hpp"
#include <vector>

#include "../GLBuffer.hpp"

class RTriangles : public RenderObj
{
public:
  RTriangles(bool hostTransfers);
  ~RTriangles();

  virtual void glRender();

  void setGLColors(std::vector<float>& VertexColor);
  void setGLPositions(std::vector<float>& VertexPos);
  void setGLNormals(std::vector<float>& VertexNormals);
  void setGLElements(std::vector<int>& Elements);

  void initOCLVertexBuffer(cl::Context& Context);
  void initOCLColorBuffer(cl::Context& Context);
  void initOCLNormBuffer(cl::Context& Context);
  void initOCLElementBuffer(cl::Context& Context);

protected:
  GLuint _colBuff;
  size_t _colBuffSize;
  cl::GLBuffer _clbuf_Colors;

  GLuint _posBuff;
  size_t _posBuffSize;
  cl::GLBuffer _clbuf_Positions;

  GLuint _normBuff;
  size_t _normBuffSize;
  cl::GLBuffer _clbuf_Normals;

  GLuint _elementBuff;
  size_t _elementBuffSize;
  cl::GLBuffer _clbuf_Elements;
};
