#include "water_simulator/renderer/entities/water_normals.h"
#include "water_simulator/renderer/algebra.h"
namespace water_simulator::renderer::entities {

void update_water_normals(std::vector<float> &vertex_normals, std::vector<float> &face_normals,
                          const std::vector<float> &heights, std::size_t resolution, const std::vector<float> &xz,
                          const std::vector<unsigned int> &indices, const std::vector<std::size_t> &count) {
  if (heights.size() != (resolution * resolution))
    throw std::invalid_argument("Invalid heights size");

  std::size_t max_face_index = (resolution - 1) * (resolution - 1) * 2;
  std::fill(face_normals.begin(), face_normals.end(), 0.0);
  for (std::size_t face_index = 0; face_index < max_face_index; ++face_index) {
    std::size_t i = indices[face_index * 3];
    std::size_t j = indices[face_index * 3 + 1];
    std::size_t k = indices[face_index * 3 + 2];
    std::array<float, 3> vi = {xz[i * 2], heights[i], xz[i * 2 + 1]};
    std::array<float, 3> vj = {xz[j * 2], heights[j], xz[j * 2 + 1]};
    std::array<float, 3> vk = {xz[k * 2], heights[k], xz[k * 2 + 1]};

    std::array<float, 3> a = {vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]};
    std::array<float, 3> b = {vk[0] - vi[0], vk[1] - vi[1], vk[2] - vi[2]};

    std::array<float, 3> cross = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    face_normals[face_index * 3] = cross[0];
    face_normals[face_index * 3 + 1] = cross[1];
    face_normals[face_index * 3 + 2] = cross[2];
  }
  std::fill(vertex_normals.begin(), vertex_normals.end(), 0.0);
  for (std::size_t face_index = 0; face_index < max_face_index; ++face_index) {
    for (std::size_t index = 0; index < 3; ++index) {
      const std::size_t vertex_index = indices[face_index * 3 + index];
      vertex_normals[vertex_index * 3] += face_normals[face_index * 3];
      vertex_normals[vertex_index * 3 + 1] += face_normals[face_index * 3 + 1];
      vertex_normals[vertex_index * 3 + 2] += face_normals[face_index * 3 + 2];
    }
  }
  for (std::size_t index = 0; index < count.size(); ++index) {
    for (std::size_t j = 0; j < 3; ++j)
      vertex_normals[index * 3 + j] /= count[index];

    const auto normal =
        normalize({vertex_normals[index * 3], vertex_normals[index * 3 + 1], vertex_normals[index * 3 + 2]});

    vertex_normals[index * 3] = normal[0];
    vertex_normals[index * 3 + 1] = normal[1];
    vertex_normals[index * 3 + 2] = normal[2];
  }
}
} // namespace water_simulator::renderer::entities
