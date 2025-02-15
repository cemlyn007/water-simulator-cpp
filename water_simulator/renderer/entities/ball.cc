#include "water_simulator/renderer/entities/ball.h"
#include "water_simulator/renderer/gl_error_macro.h"
#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/shader_context_manager.h"
#include <GL/glew.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace water_simulator::renderer::entities {

struct BallData {
  std::vector<double> vertices;
  std::vector<unsigned int> indices;
  std::vector<double> normals;
};

constexpr BallData create_ball_mesh() {
  std::vector<double> vertices;
  std::vector<unsigned int> indices;
  std::vector<double> normals;

  constexpr double radius = 0.5;
  constexpr int sectors = 36;
  constexpr int stacks = 18;

  // Generate vertices and normals
  for (int i = 0; i <= stacks; ++i) {
    double V = static_cast<double>(i) / stacks;
    double phi = V * M_PI;

    for (int j = 0; j <= sectors; ++j) {
      double U = static_cast<double>(j) / sectors;
      double theta = U * 2 * M_PI;

      double x = radius * std::sin(phi) * std::cos(theta);
      double y = radius * std::cos(phi);
      double z = radius * std::sin(phi) * std::sin(theta);

      // Add vertex coordinates
      vertices.push_back(x);
      vertices.push_back(y);
      vertices.push_back(z);

      // Calculate normal (unit vector)
      double length = std::sqrt(x * x + y * y + z * z);
      if (length == 0.0f)
        length = 1.0f;
      normals.push_back(x / length);
      normals.push_back(y / length);
      normals.push_back(z / length);
    }
  }

  // Generate indices
  for (int i = 0; i < stacks; ++i) {
    for (int j = 0; j < sectors; ++j) {
      unsigned int first = i * (sectors + 1) + j;
      unsigned int second = first + (sectors + 1);

      indices.push_back(first + 1);
      indices.push_back(second);
      indices.push_back(first);

      indices.push_back(first + 1);
      indices.push_back(second + 1);
      indices.push_back(second);
    }
  }

  return {vertices, indices, normals};
}

Ball::Ball()
    : _shader(read_file("water_simulator/renderer/shaders/simple.vs"),
              read_file("water_simulator/renderer/shaders/simple.fs")),
      _vbo(0), _vao(0), _ebo(0) {
  BallData mesh_data = create_ball_mesh();

  std::vector<float> interleaved(mesh_data.vertices.size() + mesh_data.normals.size());
  size_t dst = 0;
  for (size_t i = 0; i < mesh_data.vertices.size(); i += 3) {
    interleaved[dst + 0] = static_cast<float>(mesh_data.vertices[i]);
    interleaved[dst + 1] = static_cast<float>(mesh_data.vertices[i + 1]);
    interleaved[dst + 2] = static_cast<float>(mesh_data.vertices[i + 2]);
    interleaved[dst + 3] = static_cast<float>(mesh_data.normals[i]);
    interleaved[dst + 4] = static_cast<float>(mesh_data.normals[i + 1]);
    interleaved[dst + 5] = static_cast<float>(mesh_data.normals[i + 2]);
    dst += 6;
  }

  _vbo = init_vbo(interleaved);
  _ebo = init_ebo(mesh_data.indices);
  _vao = init_vao(_vbo, _ebo);
  glBindVertexArray(0);

  _indices = mesh_data.indices.size();
}

Ball::~Ball() {
  glBindVertexArray(0);
  if (_vbo != 0)
    glDeleteBuffers(1, &_vbo);
  if (_ebo != 0)
    glDeleteBuffers(1, &_ebo);
  if (_vao != 0)
    glDeleteVertexArrays(1, &_vao);
}

GLuint Ball::init_vbo(const std::vector<float> &vertices) {
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
  return vbo;
}

GLuint Ball::init_ebo(const std::vector<unsigned int> &indices) {
  GLuint ebo;
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
  return ebo;
}

GLuint Ball::init_vao(GLuint vbo, GLuint ebo) {
  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

  GL_CALL(glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(0)));
  glEnableVertexAttribArray(0);

  GL_CALL(
      glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(3 * sizeof(float))));
  glEnableVertexAttribArray(1);

  glBindVertexArray(0);
  return vao;
}

void Ball::set_view(const std::array<float, 16> &view) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("view", view);
  }
}

void Ball::set_view_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("viewPos", position);
  }
}

void Ball::set_projection(const std::array<float, 16> &projection) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("projection", projection);
  }
}

void Ball::set_model(const std::array<float, 16> &model) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("model", model);
  }
}

void Ball::set_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("objectColor", color);
  }
}

void Ball::set_light_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("lightPos", position);
  }
}

void Ball::set_light_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("lightColor", color);
  }
}

void Ball::draw() {
  ShaderContextManager context(_shader);
  {
    glBindVertexArray(_vao);
    glDrawElements(GL_TRIANGLES, _indices, GL_UNSIGNED_INT, nullptr);
  }
}
} // namespace water_simulator::renderer::entities