#pragma once
#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/texture.h"
#include <GL/glew.h>
#include <array>
#include <vector>

namespace water_simulator::renderer::entities {

class Water {
public:
  Water(size_t resolution, float size, float xz_offset);
  ~Water();

  void set_view(const std::array<float, 16> &view);
  void set_view_position(const std::array<float, 3> &position);
  void set_projection(const std::array<float, 16> &projection);
  void set_model(const std::array<float, 16> &model);
  void set_color(const std::array<float, 3> &color);
  void set_light_position(const std::array<float, 3> &position);
  void set_light_color(const std::array<float, 3> &color);
  void set_texture(Texture &texture);

  void set_heights(const std::vector<float> &heights);
  void set_normals(const std::vector<float> &normals);

  void draw();

private:
  size_t _resolution;

  Shader _shader;
  GLuint _xz_vbo;
  GLuint _y_vbo;
  GLuint _normal_vbo;
  GLuint _vao;
  GLuint _ebo;

  size_t _indices;

  GLuint init_vbo(const std::vector<float> &vertices);
  GLuint init_vbo(size_t bytes, bool dynamic);
  GLuint init_ebo(const std::vector<unsigned int> &indices);
  GLuint init_vao(GLuint xz_vbo, GLuint y_vbo, GLuint normal_vbo, GLuint ebo, const std::vector<float> &vertices);
};

} // namespace water_simulator::renderer::entities
