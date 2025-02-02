
#include "water_simulator/renderer/entities/light.h"
#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/shader_context_manager.h"

#include <GL/glew.h>

#include <string>

namespace water_simulator::renderer::entities {
Light::Light()
    : _shader(read_file("water_simulator/renderer/shaders/light_cube.vs"),
              read_file("water_simulator/renderer/shaders/light_cube.fs")),
      _vbo(0), _vao(0), _ebo(0) {
  std::array<float, 72> vertices = {
      // Front Face
      -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
      // Back Face
      0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
      // Top Face
      -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
      // Bottom Face
      -0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, 0.5, -0.5, -0.5, 0.5,
      // Right Face
      0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5,
      // Left Face
      -0.5, -0.5, -0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5,
      //
  };

  std::array<unsigned int, 36> indices = {
      // #Front Face
      0, 1, 2, 2, 3, 0,
      // #Back Face
      4, 5, 6, 6, 7, 4,
      // #Top Face
      8, 9, 10, 10, 11, 8,
      // #Bottom Face
      12, 13, 14, 14, 15, 12,
      // #Right Face
      16, 17, 18, 18, 19, 16,
      // #Left Face
      20, 21, 22, 22, 23, 20,
      //
  };

  _vbo = init_vbo(vertices);
  _ebo = init_ebo(indices);
  _vao = init_vao(_vbo, _ebo, vertices);

  glBindVertexArray(0);
}

Light::~Light() {
  if (_vbo != 0)
    glDeleteBuffers(1, &_vbo);
  if (_ebo != 0)
    glDeleteBuffers(1, &_ebo);
  if (_vao != 0)
    glDeleteVertexArrays(1, &_vao);
}

GLuint Light::init_vbo(const std::array<float, 72> &vertices) {
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float),
               vertices.data(), GL_STATIC_DRAW);
  return vbo;
}

GLuint Light::init_ebo(const std::array<unsigned int, 36> &indices) {
  GLuint ebo;
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
               indices.data(), GL_STATIC_DRAW);
  return ebo;
}

GLuint Light::init_vao(GLuint vbo, GLuint ebo,
                       const std::array<float, 72> &vertices) {
  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float),
                        reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);

  glBindVertexArray(0);
  return vao;
}

void Light::set_view(const std::array<float, 16> &view) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("view", view);
  }
}

void Light::set_projection(const std::array<float, 16> &projection) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("projection", projection);
  }
}

void Light::set_model(const std::array<float, 16> &model) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("model", model);
  }
}

void Light::set_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("objectColor", color);
  }
}

void Light::draw() {
  ShaderContextManager context(_shader);
  {
    glBindVertexArray(_vao);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, nullptr);
  }
}
} // namespace water_simulator::renderer::entities