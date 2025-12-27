#include "water_simulator/renderer/entities/container.h"
#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/shader_context_manager.h"

#include <GL/glew.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace water_simulator::renderer::entities {

struct MeshData {
  std::vector<float> vertices;
  std::vector<unsigned int> indices;
};

struct CubeData {
  std::vector<float> vertices;
  std::vector<float> normals;
  std::vector<unsigned int> indices;
};

CubeData cube_vertices_normals_and_indices() {
  return {// Vertices (8 corners × 3 faces each = 24 entries)
          .vertices =
              {// Front face
               0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
               // Back face
               1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f,
               // Top face
               0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
               // Bottom face
               0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
               // Right face
               1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 1.0f,
               // Left face
               0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f},

          // Normals (6 faces × 4 vertices each = 24 entries)
          .normals =
              {// Front (Z+)
               0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1,
               // Back (Z-)
               0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1,
               // Top (Y+)
               0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
               // Bottom (Y-)
               0, -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0,
               // Right (X+)
               1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
               // Left (X-)
               -1, 0, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0},

          // Indices (6 faces × 2 triangles × 3 vertices = 36 indices)
          .indices = {
              0,  1,  2,  2,  3,  0,  // Front
              4,  5,  6,  6,  7,  4,  // Back
              8,  9,  10, 10, 11, 8,  // Top
              12, 13, 14, 14, 15, 12, // Bottom
              16, 17, 18, 18, 19, 16, // Right
              20, 21, 22, 22, 23, 20  // Left
          }};
}

MeshData create_mesh(float size, float wall_thickness) {
  constexpr float height_scale = 1.25f;
  CubeData cube = cube_vertices_normals_and_indices();
  MeshData mesh;

  // Add the floor plane:
  mesh.vertices = {0, 0, 0, 0, size, 0, size, 0, 0, 0, size, 0, size, 0, size, 0, size, 0, 0, 0, size, 0, size, 0};
  mesh.indices = {2, 1, 0, 0, 3, 2};

  const std::array walls = {
      std::pair{std::array{size, height_scale, wall_thickness}, std::array{0.0f, 0.0f, size}},
      std::pair{std::array{wall_thickness, height_scale, size}, std::array{size, 0.0f, 0.0f}},

      std::pair{std::array{size, height_scale, wall_thickness}, std::array{0.0f, 0.0f, -wall_thickness}},
      std::pair{std::array{wall_thickness, height_scale, size}, std::array{-wall_thickness, 0.0f, 0.0f}},
  };
  for (const auto &[scaling, translation] : walls) {
    for (std::size_t i = 0; i < cube.vertices.size(); i += 3) {
      const float x = cube.vertices[i] * scaling[0] + translation[0];
      const float y = cube.vertices[i + 1] * scaling[1] + translation[1];
      const float z = cube.vertices[i + 2] * scaling[2] + translation[2];
      mesh.vertices.insert(mesh.vertices.end(), {x, y, z});
      mesh.vertices.insert(mesh.vertices.end(), cube.normals.begin() + i, cube.normals.begin() + i + 3);
    }
  }

  // Generate indices
  const unsigned int offset = *std::ranges::max_element(mesh.indices) + 1;
  const unsigned int max_cube_idx = *std::ranges::max_element(cube.indices);

  for (int i = 0; i < 4; ++i) {
    const unsigned int base = offset + i * (max_cube_idx + 1);
    for (const auto idx : cube.indices) {
      mesh.indices.push_back(base + idx);
    }
  }
  return mesh;
}

Container::Container(float wall_size, float wall_thickness)
    : _shader(read_file("water_simulator/renderer/shaders/simple.vs"),
              read_file("water_simulator/renderer/shaders/simple.fs")),
      _vbo(0), _vao(0), _ebo(0) {
  MeshData mesh_data = create_mesh(wall_size, wall_thickness);

  _vbo = init_vbo(mesh_data.vertices);
  _ebo = init_ebo(mesh_data.indices);
  _vao = init_vao(_vbo, _ebo, mesh_data.vertices);

  glBindVertexArray(0);
}

Container::~Container() {
  glBindVertexArray(0);
  if (_vbo != 0)
    glDeleteBuffers(1, &_vbo);
  if (_ebo != 0)
    glDeleteBuffers(1, &_ebo);
  if (_vao != 0)
    glDeleteVertexArrays(1, &_vao);
}

GLuint Container::init_vbo(const std::vector<float> &vertices) {
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
  return vbo;
}

GLuint Container::init_ebo(const std::vector<unsigned int> &indices) {
  GLuint ebo;
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
  return ebo;
}

GLuint Container::init_vao(GLuint vbo, GLuint ebo, const std::vector<float> &vertices) {
  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);

  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);

  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), reinterpret_cast<void *>(3 * sizeof(float)));
  glEnableVertexAttribArray(1);

  glBindVertexArray(0);
  return vao;
}

void Container::set_view(const std::array<float, 16> &view) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("view", view);
  }
}

void Container::set_view_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("viewPos", position);
  }
}

void Container::set_projection(const std::array<float, 16> &projection) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("projection", projection);
  }
}

void Container::set_model(const std::array<float, 16> &model) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_matrix("model", model);
  }
}

void Container::set_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("objectColor", color);
  }
}

void Container::set_light_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("lightPos", position);
  }
}

void Container::set_light_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  {
    _shader.set_uniform_vector("lightColor", color);
  }
}

void Container::draw() {
  ShaderContextManager context(_shader);
  {
    glBindVertexArray(_vao);
    glDrawElements(GL_TRIANGLES, 150, GL_UNSIGNED_INT, nullptr);
  }
}
} // namespace water_simulator::renderer::entities