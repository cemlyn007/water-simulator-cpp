#pragma once
#include <array>
#include <iostream>
#include <iterator>
namespace water_simulator::renderer {

template <std::size_t N> void print(const std::array<float, N> &arr) {
  std::copy(arr.begin(), arr.end(), std::ostream_iterator<float>(std::cout, " "));
  std::cout << "\n";
}

float radians(float degrees);

float norm(std::array<float, 3> vector);

std::array<float, 3> normalize(const std::array<float, 3> &vector);

std::array<float, 3> update_orbit_camera_position(float azimuth_radians, float elevation_radians, float radius);

std::array<float, 16> eye4d();

std::array<float, 16> multiply_matrices(const std::array<float, 16> &a, const std::array<float, 16> &b);

std::array<float, 4> multiply_matrix(const std::array<float, 16> &a, const std::array<float, 4> &b);

std::array<float, 16> transpose(const std::array<float, 16> &matrix);

std::array<float, 16> translate(const std::array<float, 16> &matrix, const std::array<float, 3> &vector);

std::array<float, 16> scale(const std::array<float, 16> &matrix, const std::array<float, 3> &vector);

std::array<float, 16> look_at(const std::array<float, 3> &eye, const std::array<float, 3> &center,
                              const std::array<float, 3> &up);

std::array<float, 16> perspective(float fov, float aspect, float near, float far);

std::array<float, 16> inverse(std::array<float, 16> matrix);

} // namespace water_simulator::renderer
