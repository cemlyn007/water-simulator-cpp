#include "water_simulator/renderer/entities/water.h"
#include "water_simulator/renderer/algebra.h"
#include "water_simulator/renderer/gl_error_macro.h"
#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/shader_context_manager.h"
#include <GL/glew.h>
#include <array>
#include <string>
#include <vector>
namespace water_simulator::renderer::entities {

struct WaterData {
  std::vector<float> vertices;
  std::vector<float> normals;
  std::vector<unsigned int> indices;
};

WaterData grid_vertices_normals_and_indices(int n_cells_x, int n_cells_z, double cell_size) {
  std::vector<float> vertices;
  std::vector<float> normals;
  std::vector<unsigned int> indices;
  if (n_cells_x <= 0 || n_cells_z <= 0 || cell_size <= 0.0f) {
    throw std::invalid_argument("Invalid grid parameters - must be positive values");
  }
  const int vertices_x = n_cells_x;
  const int vertices_z = n_cells_z;
  const int total_vertices = vertices_x * vertices_z;
  vertices.reserve(total_vertices);
  normals.reserve(total_vertices);
  indices.reserve(n_cells_x * n_cells_z * 6);
  for (size_t x = 0; x < vertices_x; ++x) {
    for (size_t z = 0; z < vertices_z; ++z) {
      vertices.push_back(x * cell_size);
      vertices.push_back(z * cell_size);
      normals.push_back(0.0f);
      normals.push_back(1.0f);
      normals.push_back(0.0f);
    }
  }

  for (int x = 0; x < n_cells_x - 1; ++x) {
    for (int z = 0; z < n_cells_z - 1; ++z) {
      unsigned int top_left = x * vertices_z + z;
      unsigned int top_right = top_left + 1;
      unsigned int bottom_left = (x + 1) * vertices_z + z;
      unsigned int bottom_right = bottom_left + 1;
      indices.push_back(top_left);
      indices.push_back(top_right);
      indices.push_back(bottom_right);
      indices.push_back(bottom_right);
      indices.push_back(bottom_left);
      indices.push_back(top_left);
    }
  }
  return {vertices, normals, indices};
}

Water::Water(size_t resolution, float length, float xz_offset)
    : _resolution(resolution), _shader(read_file("water_simulator/renderer/shaders/basic_lighting.vs"),
                                       read_file("water_simulator/renderer/shaders/basic_lighting.fs")),
      _xz_vbo(0), _y_vbo(0), _normal_vbo(0), _vao(0), _ebo(0), _vertex_normals(3 * _resolution * _resolution),
      _face_normals((_resolution - 1) * (_resolution - 1) * 2 * 3), _count(_resolution * _resolution, 0) {
  WaterData mesh_data = grid_vertices_normals_and_indices(_resolution, _resolution, length / (_resolution));
  std::transform(mesh_data.vertices.begin(), mesh_data.vertices.end(), mesh_data.vertices.begin(),
                 [&](auto &value) { return value + xz_offset; });
  _xz = mesh_data.vertices;
  _xz_vbo = init_vbo(mesh_data.vertices);
  _y_vbo = init_vbo(_resolution * _resolution * sizeof(float), true);
  _normal_vbo = init_vbo(3 * _resolution * _resolution * sizeof(float), true);
  _ebo = init_ebo(mesh_data.indices);
  _vao = init_vao(_xz_vbo, _y_vbo, _normal_vbo, _ebo, mesh_data.vertices);
  _indices = mesh_data.indices;
  glBindVertexArray(0);

  size_t max_face_index = (_resolution - 1) * (_resolution - 1) * 2;
  for (size_t face_index = 0; face_index < max_face_index; ++face_index) {
    for (size_t index = 0; index < 3; ++index) {
      const size_t vertex_index = _indices[face_index * 3 + index];
      ++_count[vertex_index];
    }
  }
}

Water::~Water() {
  if (_xz_vbo != 0)
    glDeleteBuffers(1, &_xz_vbo);
  if (_y_vbo != 0)
    glDeleteBuffers(1, &_y_vbo);
  if (_normal_vbo != 0)
    glDeleteBuffers(1, &_normal_vbo);
  if (_ebo != 0)
    glDeleteBuffers(1, &_ebo);
  if (_vao != 0)
    glDeleteVertexArrays(1, &_vao);
}

GLuint Water::init_vbo(const std::vector<float> &vertices) {
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
  return vbo;
}

GLuint Water::init_vbo(size_t bytes, bool dynamic) {
  GLuint vbo;
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, bytes, nullptr, dynamic ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
  return vbo;
}

GLuint Water::init_ebo(const std::vector<unsigned int> &indices) {
  GLuint ebo;
  glGenBuffers(1, &ebo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);
  return ebo;
}

GLuint Water::init_vao(GLuint xz_vbo, GLuint y_vbo, GLuint normal_vbo, GLuint ebo, const std::vector<float> &vertices) {
  GLuint vao;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
  glBindBuffer(GL_ARRAY_BUFFER, xz_vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, normal_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(1);
  glBindBuffer(GL_ARRAY_BUFFER, y_vbo);
  glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(float), reinterpret_cast<void *>(0));
  glEnableVertexAttribArray(2);
  glBindVertexArray(0);
  return vao;
}

void Water::set_view(const std::array<float, 16> &view) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_matrix("view", view);
}

void Water::set_view_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_vector("viewPos", position);
}

void Water::set_projection(const std::array<float, 16> &projection) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_matrix("projection", projection);
}

void Water::set_model(const std::array<float, 16> &model) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_matrix("model", model);
}

void Water::set_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_vector("objectColor", color);
}

void Water::set_light_position(const std::array<float, 3> &position) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_vector("lightPos", position);
}

void Water::set_light_color(const std::array<float, 3> &color) {
  ShaderContextManager context(_shader);
  _shader.set_uniform_vector("lightColor", color);
}

void Water::set_texture(Texture &texture) {
  texture.use();
  {
    ShaderContextManager context(_shader);
    _shader.set_uniform("background", 0);
  }
}

void Water::set_heights(const std::vector<float> &heights) {
  if (heights.size() != (_resolution * _resolution))
    throw std::invalid_argument("Invalid heights size");
  glBindBuffer(GL_ARRAY_BUFFER, _y_vbo);
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, heights.size() * sizeof(float), heights.data(), GL_DYNAMIC_DRAW));
  update_normals(heights);
  set_normals(_vertex_normals);
}

void Water::set_normals(const std::vector<float> &normals) {
  if (normals.size() != (3 * _resolution * _resolution))
    throw std::invalid_argument("Invalid normals size");
  glBindBuffer(GL_ARRAY_BUFFER, _normal_vbo);
  GL_CALL(glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_DYNAMIC_DRAW));
}

void Water::draw() {
  ShaderContextManager context(_shader);
  glBindVertexArray(_vao);
  glDrawElements(GL_TRIANGLES, _indices.size(), GL_UNSIGNED_INT, nullptr);
}

void Water::update_normals(const std::vector<float> &heights) {
  if (heights.size() != (_resolution * _resolution))
    throw std::invalid_argument("Invalid heights size");

  size_t max_face_index = (_resolution - 1) * (_resolution - 1) * 2;
  std::fill(_face_normals.begin(), _face_normals.end(), 0.0);
  for (size_t face_index = 0; face_index < max_face_index; ++face_index) {
    size_t i = _indices[face_index * 3];
    size_t j = _indices[face_index * 3 + 1];
    size_t k = _indices[face_index * 3 + 2];
    std::array<float, 3> vi = {_xz[i * 2], heights[i], _xz[i * 2 + 1]};
    std::array<float, 3> vj = {_xz[j * 2], heights[j], _xz[j * 2 + 1]};
    std::array<float, 3> vk = {_xz[k * 2], heights[k], _xz[k * 2 + 1]};

    std::array<float, 3> a = {vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]};
    std::array<float, 3> b = {vk[0] - vi[0], vk[1] - vi[1], vk[2] - vi[2]};

    std::array<float, 3> cross = {a[1] * b[2] - a[2] * b[1] + a[1] * b[2] - a[2] * b[1],
                                  a[2] * b[0] - a[0] * b[2] + a[2] * b[0] - a[0] * b[2],
                                  a[0] * b[1] - a[1] * b[0] + a[0] * b[1] - a[1] * b[0]};
    _face_normals[face_index * 3] = cross[0];
    _face_normals[face_index * 3 + 1] = cross[1];
    _face_normals[face_index * 3 + 2] = cross[2];
  }
  std::fill(_vertex_normals.begin(), _vertex_normals.end(), 0.0);
  for (size_t face_index = 0; face_index < max_face_index; ++face_index) {
    for (size_t index = 0; index < 3; ++index) {
      const size_t vertex_index = _indices[face_index * 3 + index];
      _vertex_normals[vertex_index * 3] += _face_normals[face_index * 3];
      _vertex_normals[vertex_index * 3 + 1] += _face_normals[face_index * 3 + 1];
      _vertex_normals[vertex_index * 3 + 2] += _face_normals[face_index * 3 + 2];
    }
  }
  for (size_t index = 0; index < _count.size(); ++index) {
    for (size_t j = 0; j < 3; ++j)
      _vertex_normals[index * 3 + j] /= _count[index];

    const auto normal =
        normalize({_vertex_normals[index * 3], _vertex_normals[index * 3 + 1], _vertex_normals[index * 3 + 2]});

    _vertex_normals[index * 3] = normal[0];
    _vertex_normals[index * 3 + 1] = normal[1];
    _vertex_normals[index * 3 + 2] = normal[2];
  }
}

} // namespace water_simulator::renderer::entities