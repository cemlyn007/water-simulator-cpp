#pragma once
#include "water_simulator/renderer/shader.h"
#include <GL/glew.h>
#include <array>

namespace water_simulator::renderer::entities {

class Light {
public:
  Light();
  ~Light();

  void set_view(const std::array<float, 16> &view);
  void set_projection(const std::array<float, 16> &projection);
  void set_model(const std::array<float, 16> &model);
  void set_color(const std::array<float, 3> &color);

  void draw();

private:
  Shader _shader;
  GLuint _vbo;
  GLuint _vao;
  GLuint _ebo;

  GLuint init_vbo(const std::array<float, 72> &vertices);
  GLuint init_ebo(const std::array<unsigned int, 36> &indices);
  GLuint init_vao(GLuint vbo, GLuint ebo, const std::array<float, 72> &vertices);
};

} // namespace water_simulator::renderer::entities
