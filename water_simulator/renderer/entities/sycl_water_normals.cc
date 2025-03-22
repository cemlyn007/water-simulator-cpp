#include "water_simulator/renderer/entities/water_normals.h"
#include <cmath>
#include <sycl/sycl.hpp>
namespace water_simulator::renderer::entities {

static sycl::queue queue;

void update_water_normals(std::vector<float> &vertex_normals, std::vector<float> &face_normals,
                          const std::vector<float> &heights, size_t resolution, const std::vector<float> &xz,
                          const std::vector<unsigned int> &indices, const std::vector<size_t> &count) {
  if (heights.size() != (resolution * resolution))
    throw std::invalid_argument("Invalid heights size");

  const size_t null_face_index = (indices.size() / 3) - 1;
  const size_t vertex_size = vertex_normals.size() / 3;
  const size_t max_face_index = (resolution - 1) * (resolution - 1) * 2;
  constexpr size_t max_faces_per_vertex = 6;

  static std::vector<float> vertex_face_indices(vertex_normals.size() * max_faces_per_vertex, null_face_index);
  static std::vector<size_t> vertex_face_counts(vertex_normals.size(), 0);

  if (face_normals.size() != (resolution - 1) * (resolution - 1) * 2 * 3 + 3) {
    // The +1 is cheeky, it's going to act as the null face index
    face_normals.push_back(0.0);
    face_normals.push_back(0.0);
    face_normals.push_back(0.0);
    for (size_t face_index = 0; face_index < indices.size() / 3; ++face_index) {
      for (size_t offset = 0; offset < 3; ++offset) {
        const size_t vertex_index = indices[face_index * 3 + offset];
        const size_t space = vertex_face_counts[vertex_index];
        vertex_face_indices[vertex_index * max_faces_per_vertex + space] = face_index;
        vertex_face_counts[vertex_index]++;
      }
    }
  }

  sycl::buffer<float, 1> d_face_normals(face_normals.data(), sycl::range<1>(face_normals.size()));
  sycl::buffer<float, 1> d_heights(heights.data(), sycl::range<1>(heights.size()));
  sycl::buffer<float, 1> d_xz(xz.data(), sycl::range<1>(xz.size()));
  sycl::buffer<unsigned int, 1> d_indices(indices.data(), sycl::range<1>(indices.size()));

  sycl::buffer<float, 1> d_vertex_normals(vertex_normals.data(), sycl::range<1>(vertex_normals.size()));
  sycl::buffer<float, 1> d_vertex_face_indices(vertex_face_indices.data(), sycl::range<1>(vertex_face_indices.size()));
  sycl::buffer<size_t, 1> d_vertex_face_counts(vertex_face_counts.data(), sycl::range<1>(vertex_face_counts.size()));

  auto event = queue.submit([&](sycl::handler &handler) {
    // Left as write because the last face normal needs to stay as all zeros in a later step.
    auto face_normals_acc = d_face_normals.get_access<sycl::access::mode::write>(handler);
    auto heights_acc = d_heights.get_access<sycl::access::mode::read>(handler);
    auto xz_acc = d_xz.get_access<sycl::access::mode::read>(handler);
    auto indices_acc = d_indices.get_access<sycl::access::mode::read>(handler);
    handler.parallel_for(sycl::range<1>(max_face_index), [=](sycl::id<1> index) {
      const size_t face_index = index;
      const size_t i = indices_acc[face_index * 3];
      const size_t j = indices_acc[face_index * 3 + 1];
      const size_t k = indices_acc[face_index * 3 + 2];
      std::array<float, 3> vi = {xz_acc[i * 2], heights_acc[i], xz_acc[i * 2 + 1]};
      std::array<float, 3> vj = {xz_acc[j * 2], heights_acc[j], xz_acc[j * 2 + 1]};
      std::array<float, 3> vk = {xz_acc[k * 2], heights_acc[k], xz_acc[k * 2 + 1]};

      std::array<float, 3> a = {vj[0] - vi[0], vj[1] - vi[1], vj[2] - vi[2]};
      std::array<float, 3> b = {vk[0] - vi[0], vk[1] - vi[1], vk[2] - vi[2]};

      std::array<float, 3> cross = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
      face_normals_acc[face_index * 3] = cross[0];
      face_normals_acc[face_index * 3 + 1] = cross[1];
      face_normals_acc[face_index * 3 + 2] = cross[2];
    });
  });
  queue.submit([&](sycl::handler &handler) {
    handler.depends_on(event);
    auto vertex_normals_acc = d_vertex_normals.get_access<sycl::access::mode::discard_write>(handler);
    auto face_normals_acc = d_face_normals.get_access<sycl::access::mode::read>(handler);
    auto vertex_face_indices_acc = d_vertex_face_indices.get_access<sycl::access::mode::read>(handler);
    auto vertex_face_counts_acc = d_vertex_face_counts.get_access<sycl::access::mode::read>(handler);
    handler.parallel_for(sycl::range<1>(vertex_size), [=](sycl::id<1> index) {
      auto vertex_index = index;
      float x = 0.0;
      float y = 0.0;
      float z = 0.0;
      const size_t faces = vertex_face_counts_acc[vertex_index];
      for (size_t offset = 0; offset < max_faces_per_vertex; ++offset) {
        const size_t face_index = vertex_face_indices_acc[vertex_index * max_faces_per_vertex + offset];
        x += face_normals_acc[face_index * 3];
        y += face_normals_acc[face_index * 3 + 1];
        z += face_normals_acc[face_index * 3 + 2];
      }
      x /= faces;
      y /= faces;
      z /= faces;
      float norm_squared = x * x + y * y + z * z;
      if (norm_squared == 0.0) {
        x = 0.0;
        y = 1.0;
        z = 0.0;
      } else {
        float norm = std::sqrt(norm_squared);
        x /= norm;
        y /= norm;
        z /= norm;
      }
      vertex_normals_acc[vertex_index * 3] = x;
      vertex_normals_acc[vertex_index * 3 + 1] = y;
      vertex_normals_acc[vertex_index * 3 + 2] = z;
    });
  });
  queue.wait();
}
} // namespace water_simulator::renderer::entities
