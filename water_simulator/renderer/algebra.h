#pragma once
#include <array>

namespace water_simulator::renderer {

void print(const std::array<float, 3> &vector);

void print(const std::array<float, 16> &matrix);

void print(const std::array<std::array<float, 4>, 4> &matrix);

float radians(float degrees);

float norm(std::array<float, 3> vector);

std::array<float, 3> update_orbit_camera_position(float azimuth_radians, float elevation_radians, float radius);

std::array<float, 16> eye4d();

std::array<float, 16> multiply_matrices(const std::array<float, 16> &a, const std::array<float, 16> &b);

std::array<float, 16> translate(const std::array<float, 16> &matrix, const std::array<float, 3> &vector);

std::array<float, 16> scale(const std::array<float, 16> &matrix, const std::array<float, 3> &vector);

std::array<float, 16> look_at(const std::array<float, 3> &eye, const std::array<float, 3> &center,
                              const std::array<float, 3> &up);

std::array<float, 16> perspective(float fov, float aspect, float near, float far);

} // namespace water_simulator::renderer
