#pragma once

#include <GL/glew.h>

#include "water_simulator/renderer/texture.h"

namespace water_simulator::renderer {
class Camera {
public:
  GLuint rendered_texture;

  Camera(int width, int height);
  ~Camera();

  void resize(int width, int height);
  void bind();
  void unbind();

  Texture texture();

private:
  int _width;
  int _height;
  GLuint _framebuffer;
  GLuint _depth_render_buffer;
};
} // namespace water_simulator::renderer
