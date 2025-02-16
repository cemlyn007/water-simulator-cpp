#include "water_simulator/engine/engine.h"
#include <cmath>
#include <mdspan>

namespace water_simulator::engine {

float sphere_mass(float radius, float density) { return 4.0 / 3.0 * M_PI * std::pow(radius, 3) * density; }

const std::vector<float> sphere_body_heights(const std::vector<float> &sphere_centers,
                                             const std::vector<float> &sphere_radii,
                                             const std::vector<float> &water_xzs,
                                             const std::vector<float> &water_heights) {
  const size_t n_spheres = sphere_centers.size() / 3;
  const size_t n_water_points = water_xzs.size() / 2;
  std::vector<float> body_heights(n_spheres * n_water_points, 0.0);
  for (size_t sphere_index = 0; sphere_index < n_spheres; ++sphere_index) {
    for (size_t water_index = 0; water_index < n_water_points; ++water_index) {
      const float x = water_xzs[2 * water_index];
      const float z = water_xzs[2 * water_index + 1];
      const float y = water_heights[water_index];
      const float sphere_x = sphere_centers[3 * sphere_index];
      const float sphere_z = sphere_centers[3 * sphere_index + 2];
      const float sphere_radius = sphere_radii[sphere_index];
      const float distance = std::sqrt(std::pow(x - sphere_x, 2) + std::pow(z - sphere_z, 2));
      if (distance < sphere_radius) {
        const float half_body_height = std::sqrt(std::pow(sphere_radius, 2) - std::pow(distance, 2));
        // If the sphere center - the half body height is above the water, then
        // the sphere is above the water.
        const float min_body_height = std::max(sphere_centers[3 * sphere_index + 1] - half_body_height, 0.0f);
        const float max_body_height = std::min(sphere_centers[3 * sphere_index + 1] + half_body_height, y);
        const float body_height = max_body_height - min_body_height;
        if (body_height > 0) {
          body_heights[sphere_index * n_water_points + water_index] = body_height;
        }
      }
    }
  }
  return body_heights;
}

const std::vector<float> cross_correlation(const std::span<float> &input, const std::vector<float> &kernel,
                                           const size_t input_n, const size_t input_m, const size_t kernel_n,
                                           const size_t kernel_m) {
  if (input.size() != input_n * input_m)
    throw std::runtime_error("Invalid input size");
  if (kernel.size() != kernel_n * kernel_m)
    throw std::runtime_error("Invalid kernel size");
  std::vector<float> output(input.size(), 0);
  for (size_t i = 0; i < input_n; ++i) {
    for (size_t j = 0; j < input_m; ++j) {
      const size_t update_index = i * input_m + j;
      for (size_t ki = 0; ki < kernel_n; ++ki) {
        size_t get_i = i + ki - kernel_n / 2;
        if (i + ki < kernel_n / 2) {
          get_i = 0;
        } else if (get_i >= input_n) {
          get_i = input_n - 1;
        }
        for (size_t kj = 0; kj < kernel_m; ++kj) {
          size_t get_j = j + kj - kernel_m / 2;
          if (j + kj < kernel_m / 2) {
            get_j = 0;
          } else if (get_j >= input_m) {
            get_j = input_m - 1;
          }
          output[update_index] += input[get_i * input_m + get_j] * kernel[ki * kernel_m + kj];
        }
      }
    }
  }
  return output;
}

std::vector<float> smooth_sphere_body_heights(std::vector<float> body_heights, const size_t n_spheres, const size_t n,
                                              const size_t m) {
  // Smoothen the body height field.
  auto span_body_heights = std::span(body_heights);
  for (size_t index = 0; index < 2; ++index) {
    for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
      auto sphere_body_heights = span_body_heights.subspan(sphere * n * m, n * m);
      const auto new_sphere_body_heights = cross_correlation(
          sphere_body_heights,
          {0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.11111111, 0.11111111},
          n, m, 3, 3);
      std::copy(new_sphere_body_heights.begin(), new_sphere_body_heights.end(), sphere_body_heights.begin());
    }
  }
  return body_heights;
}

State apply_container_collisions(State state) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    size_t sphere_y_index = 3 * sphere + 1;
    if (state._sphere_centers[sphere_y_index] - state._sphere_radii[sphere] <= 0.0) {
      state._sphere_centers[sphere_y_index] = state._sphere_radii[sphere];
      state._sphere_velocities[sphere_y_index] = 0.0;
      // TODO: I should also apply an opposite force to the sphere.
    }
  }
  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    if (state._water_heights[index] < 0.0) {
      state._water_heights[index] = 0.0;
      state._water_velocities[index] = 0.0;
    }
  }
  return state;
};

constexpr float GRAVITY = -9.81;

State apply_sphere_water_interaction(State state, const std::vector<float> &sphere_body_heights,
                                     const std::vector<float> &sphere_masses) {

  const size_t n_spheres = state._sphere_centers.size() / 3;
  const float wave_speed = std::min(state._wave_speed, 0.5f * state._spacing / state._time_delta);
  const float c = std::pow(wave_speed / state._spacing, 2);

  std::vector<float> body_heights(state._n * state._m, 0.0);
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        body_heights[i * state._m + j] += sphere_body_heights[sphere * state._m + j];
      }
    }
  }

  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_heights[index] += 0.05 * (body_heights[index] - state._body_heights[index]);
  }
  state._body_heights = body_heights;

  const float KERNEL_SUM = 4;
  auto sums = cross_correlation(state._water_heights,
                                {
                                    0.0,
                                    1.0,
                                    0.0,
                                    1.0,
                                    0.0,
                                    1.0,
                                    0.0,
                                    1.0,
                                    0.0,
                                },
                                state._n, state._m, 3, 3);

  // TODO: Export
  const float VELOCITY_DAMPING = 0.3;
  auto velocity_damping = std::max(0.0, 1.0 - VELOCITY_DAMPING * state._time_delta);
  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_velocities[index] += state._time_delta * c * (sums[index] - KERNEL_SUM * state._water_heights[index]);
    state._water_velocities[index] *= velocity_damping;
  }

  const float POSITIONAL_DAMPING = 1.0f;
  auto positional_damping = std::min(POSITIONAL_DAMPING * state._time_delta, 1.0f);
  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_heights[index] += (sums[index] / KERNEL_SUM - state._water_heights[index]) * positional_damping;
    state._water_heights[index] += state._time_delta * state._water_velocities[index];
  }

  // Sphere
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    float sphere_body_height = 0;
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        sphere_body_height += sphere_body_heights[i * state._m + j];
      }
    }
    const float force = -sphere_body_height * std::pow(state._spacing, 2) * GRAVITY;
    const float acceleration = force / sphere_masses[sphere] + GRAVITY;
    state._sphere_velocities[3 * sphere + 1] += state._time_delta * acceleration;
    if (sphere_body_height > 0)
      state._sphere_velocities[3 * sphere + 1] *= 0.999;
    state._sphere_centers[3 * sphere + 1] += state._time_delta * state._sphere_velocities[3 * sphere + 1];
  }

  return state;
}

State step(const State &state) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  if (n_spheres != state._sphere_velocities.size() / 3 || n_spheres != state._sphere_radii.size() ||
      n_spheres != state._sphere_densities.size()) {
    throw std::runtime_error("Invalid state: sphere vectors have different sizes");
  }

  std::vector<float> sphere_masses(n_spheres);
  for (size_t i = 0; i < n_spheres; ++i) {
    sphere_masses[i] = sphere_mass(state._sphere_radii[i], state._sphere_densities[i]);
  }

  State new_state = apply_sphere_water_interaction(
      state,
      smooth_sphere_body_heights(
          sphere_body_heights(state._sphere_centers, state._sphere_radii, state._water_xzs, state._water_heights),
          n_spheres, state._n, state._m),
      sphere_masses);
  new_state = apply_container_collisions(new_state);
  return new_state;
}

} // namespace water_simulator::engine