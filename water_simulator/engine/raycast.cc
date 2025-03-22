#include "water_simulator/engine/engine.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace water_simulator::engine {

std::optional<std::pair<size_t, float>> raycast(const std::vector<float> &sphere_centers,
                                                const std::vector<float> &sphere_radii,
                                                const std::array<float, 3> &ray_start,
                                                const std::array<float, 3> &ray_direction) {
  assert(sphere_centers.size() == sphere_radii.size() * 3);
  const size_t n_spheres = sphere_radii.size();
  std::optional<std::pair<size_t, float>> closest_sphere;
  float closest_distance_squared = std::numeric_limits<float>::max();

  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    auto sphere_position = std::span(sphere_centers).subspan(3 * sphere, 3);
    float sphere_radius = sphere_radii[sphere];
    float a = (ray_direction[0] * ray_direction[0] + ray_direction[1] * ray_direction[1] +
               ray_direction[2] * ray_direction[2]);
    float b = 2 * (ray_direction[0] * (ray_start[0] - sphere_position[0]) +
                   ray_direction[1] * (ray_start[1] - sphere_position[1]) +
                   ray_direction[2] * (ray_start[2] - sphere_position[2]));
    float c = sphere_position[0] * sphere_position[0] + sphere_position[1] * sphere_position[1] +
              sphere_position[2] * sphere_position[2] + ray_start[0] * ray_start[0] + ray_start[1] * ray_start[1] +
              ray_start[2] * ray_start[2] -
              2 * (sphere_position[0] * ray_start[0] + sphere_position[1] * ray_start[1] +
                   sphere_position[2] * ray_start[2]) -
              sphere_radius * sphere_radius;

    float discriminant = b * b - 4 * a * c;
    if (discriminant >= 0) {
      float scale_factor = (-b + std::sqrt(discriminant)) / (2 * a);
      float distance_squared = ((scale_factor * ray_direction[0]) * (scale_factor * ray_direction[0]) +
                                (scale_factor * ray_direction[1]) * (scale_factor * ray_direction[1]) +
                                (scale_factor * ray_direction[2]) * (scale_factor * ray_direction[2]));
      if (distance_squared < closest_distance_squared) {

        closest_distance_squared = distance_squared;
        closest_sphere = {sphere, distance_squared};
      }
    }
  }
  return closest_sphere;
}

// If the ray intersects the ground, then this returns the squared distance.
std::optional<float> raycast(const std::array<float, 3> &ray_start, const std::array<float, 3> &ray_direction) {
  const float y_direction = ray_direction[1];
  if (y_direction == 0)
    return {};
  const float lambda = -ray_start[1] / (y_direction);
  if (lambda < 0)
    return {};
  return lambda * lambda *
         (ray_direction[0] * ray_direction[0] + ray_direction[1] * ray_direction[1] +
          ray_direction[2] * ray_direction[2]);
}

std::optional<float> raycast_xz_squared(const std::array<float, 3> &ray_start,
                                        const std::array<float, 3> &ray_direction, float face_y,
                                        const std::array<float, 2> &face_start, const std::array<float, 2> &face_end) {

  const float y_direction = ray_direction[1];
  if (y_direction == 0)
    return {};
  const float lambda = (face_y - ray_start[1]) / (y_direction);
  if (lambda < 0)
    return {};
  const std::array<float, 3> intersection = {
      ray_start[0] + lambda * ray_direction[0],
      ray_start[1] + lambda * ray_direction[1],
      ray_start[2] + lambda * ray_direction[2],
  };
  const bool inbounds = ((face_start[0] < intersection[0] && intersection[0] < face_end[0]) &&
                         (face_start[1] < intersection[2] && intersection[2] < face_end[1]));
  if (!inbounds)
    return {};
  return lambda * lambda *
         (ray_direction[0] * ray_direction[0] + ray_direction[1] * ray_direction[1] +
          ray_direction[2] * ray_direction[2]);
}

std::optional<float> raycast_xy_squared(const std::array<float, 3> &ray_start,
                                        const std::array<float, 3> &ray_direction, float face_z,
                                        const std::array<float, 2> &face_start, const std::array<float, 2> &face_end) {

  const float z_direction = ray_direction[2];
  if (z_direction == 0)
    return {};
  const float lambda = (face_z - ray_start[2]) / (z_direction);
  if (lambda < 0)
    return {};
  const std::array<float, 3> intersection = {
      ray_start[0] + lambda * ray_direction[0],
      ray_start[1] + lambda * ray_direction[1],
      ray_start[2] + lambda * ray_direction[2],
  };
  const bool inbounds = ((face_start[0] < intersection[1] && intersection[1] < face_end[0]) &&
                         (face_start[1] < intersection[2] && intersection[2] < face_end[1]));
  if (!inbounds)
    return {};
  return lambda * lambda *
         (ray_direction[0] * ray_direction[0] + ray_direction[1] * ray_direction[1] +
          ray_direction[2] * ray_direction[2]);
}

std::optional<float> raycast_yz_squared(const std::array<float, 3> &ray_start,
                                        const std::array<float, 3> &ray_direction, float face_x,
                                        const std::array<float, 2> &face_start, const std::array<float, 2> &face_end) {

  const float x_direction = ray_direction[0];
  if (x_direction == 0)
    return {};
  const float lambda = (face_x - ray_start[0]) / (x_direction);
  if (lambda < 0)
    return {};
  const std::array<float, 3> intersection = {
      ray_start[0] + lambda * ray_direction[0],
      ray_start[1] + lambda * ray_direction[1],
      ray_start[2] + lambda * ray_direction[2],
  };
  const bool inbounds = ((face_start[0] < intersection[1] && intersection[1] < face_end[0]) &&
                         (face_start[1] < intersection[2] && intersection[2] < face_end[1]));
  if (!inbounds)
    return {};
  return lambda * lambda *
         (ray_direction[0] * ray_direction[0] + ray_direction[1] * ray_direction[1] +
          ray_direction[2] * ray_direction[2]);
}

std::optional<std::pair<size_t, float>> raycast(const std::vector<float> &sphere_centers,
                                                const std::vector<float> &sphere_radii,
                                                const std::array<float, 3> &ray_start,
                                                const std::array<float, 3> &ray_direction, const float wall_size,
                                                const float wall_thickness, const float wall_height) {
  assert(sphere_centers.size() == sphere_radii.size() * 3);
  auto closest_sphere_result = raycast(sphere_centers, sphere_radii, ray_start, ray_direction);
  if (!closest_sphere_result.has_value()) {
    return {};
  }
  const float sphere_distance_squared = closest_sphere_result.value().second;

  const auto floor_distance_squared = raycast(ray_start, ray_direction);
  if (floor_distance_squared.has_value() && floor_distance_squared.value() < sphere_distance_squared)
    return {};

  const float min_point = -wall_size / 2.0f;
  const float max_point = wall_size / 2.0f;

  // North Wall Interior Face
  auto distance_squared =
      raycast_xy_squared(ray_start, ray_direction, max_point, {min_point, 0.0}, {max_point, wall_height});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // North Wall Exterior Face
  distance_squared = raycast_xy_squared(ray_start, ray_direction, max_point + wall_thickness, {min_point, 0.0},
                                        {max_point, wall_height});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // North Wall Top Face
  distance_squared = raycast_xz_squared(ray_start, ray_direction, wall_height, {min_point, max_point},
                                        {max_point, max_point + wall_thickness});
  if (distance_squared.has_value())
    std::cout << distance_squared.value() << std::endl;
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // South Wall Interior Face
  distance_squared =
      raycast_xy_squared(ray_start, ray_direction, min_point, {min_point, 0.0}, {max_point, wall_height});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // South Wall Exterior Face
  distance_squared = raycast_xy_squared(ray_start, ray_direction, min_point - wall_thickness, {min_point, 0.0},
                                        {max_point, wall_height});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // South Wall Top Face
  distance_squared = raycast_xz_squared(ray_start, ray_direction, wall_height, {min_point, min_point - wall_thickness},
                                        {max_point, min_point});
  if (distance_squared.has_value())
    std::cout << distance_squared.value() << std::endl;
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // East Wall Interior Face
  distance_squared =
      raycast_yz_squared(ray_start, ray_direction, max_point, {0.0, min_point}, {wall_height, max_point});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // East Wall Exterior Face
  distance_squared = raycast_yz_squared(ray_start, ray_direction, max_point + wall_thickness, {0.0, min_point},
                                        {wall_height, max_point});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // East Wall Top Face
  distance_squared = raycast_xz_squared(ray_start, ray_direction, wall_height, {max_point, min_point},
                                        {max_point, max_point + wall_thickness});
  if (distance_squared.has_value())
    std::cout << distance_squared.value() << std::endl;
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // West Wall Interior Face
  distance_squared =
      raycast_yz_squared(ray_start, ray_direction, min_point, {0.0, min_point}, {wall_height, max_point});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // West Wall Exterior Face
  distance_squared = raycast_yz_squared(ray_start, ray_direction, min_point - wall_thickness, {0.0, min_point},
                                        {wall_height, max_point});
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  // West Wall Top Face
  distance_squared = raycast_xz_squared(ray_start, ray_direction, wall_height, {min_point - wall_thickness, min_point},
                                        {min_point, max_point});
  if (distance_squared.has_value())
    std::cout << distance_squared.value() << std::endl;
  if (distance_squared.has_value() && distance_squared.value() < sphere_distance_squared)
    return {};

  return closest_sphere_result;
}

} // namespace water_simulator::engine