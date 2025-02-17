#include "water_simulator/engine/engine.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mdspan>
#include <vector>

namespace water_simulator::engine {

void saveVectorToCSV(const std::vector<float> &data, const std::string &filename) {
  std::ofstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    return;
  }

  for (size_t i = 0; i < data.size(); ++i) {
    file << data[i];
    if (i < data.size() - 1) {
      file << ","; // Separate values with commas
    }
  }
  file << "\n"; // End with a newline
  file.close();

  std::cout << "Data successfully written to " << filename << std::endl;
}

float sphere_mass(float radius, float density) { return 4.0 / 3.0 * M_PI * std::pow(radius, 3) * density; }

const std::vector<float> sphere_body_heights(const std::vector<float> &sphere_centers,
                                             const std::vector<float> &sphere_radii,
                                             const std::vector<float> &water_xzs,
                                             const std::vector<float> &water_heights) {
  assert(sphere_centers.size() % 3 == 0);
  assert(water_xzs.size() % 2 == 0);
  assert(sphere_radii.size() * 3 == sphere_centers.size());

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
        int get_i = i + ki - kernel_n / 2;
        if (i + ki < kernel_n / 2) {
          get_i = 0;
        } else if (get_i >= input_n) {
          get_i = input_n - 1;
        }
        for (size_t kj = 0; kj < kernel_m; ++kj) {
          int get_j = j + kj - kernel_m / 2;
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

State apply_container_collisions(State state, float restitution) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    size_t sphere_y_index = 3 * sphere + 1;
    if (state._sphere_centers[sphere_y_index] <= state._sphere_radii[sphere]) {
      state._sphere_centers[sphere_y_index] = state._sphere_radii[sphere];
      state._sphere_velocities[sphere_y_index] = -state._sphere_velocities[sphere_y_index] * restitution;
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

  assert(state._water_heights.size() == state._n * state._m);
  assert(sphere_body_heights.size() == n_spheres * state._n * state._m);
  assert(state._water_velocities.size() == state._water_heights.size());

  std::vector<float> body_heights(state._n * state._m, 0.0);
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        body_heights[i * state._m + j] += sphere_body_heights[sphere * state._n * state._m + i * state._m + j];
      }
    }
  }

  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_heights[index] += 0.05 * (body_heights[index] - state._body_heights[index]);
  }
  state._body_heights = body_heights;

  const float KERNEL_SUM = 4;
  // state._water_heights = std::vector<float>(state._n * state._m, 1.0);
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

  // saveVectorToCSV(sums, "/home/cemlyn/Development/water-simulator-cpp/sums.csv");
  // throw std::runtime_error("Not implemented");

  const float wave_speed = std::min(state._wave_speed, 0.5f * state._spacing / state._time_delta);
  const float c = std::pow(wave_speed / state._spacing, 2);

  // TODO: Export
  const float POSITIONAL_DAMPING = 1.0f;
  auto positional_damping = std::min(POSITIONAL_DAMPING * state._time_delta, 1.0f);
  const float VELOCITY_DAMPING = 0.3;
  auto velocity_damping = std::max(0.0, 1.0 - VELOCITY_DAMPING * state._time_delta);
  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_velocities[index] += state._time_delta * c * (sums[index] - KERNEL_SUM * state._water_heights[index]);
    state._water_heights[index] += (sums[index] / KERNEL_SUM - state._water_heights[index]) * positional_damping;
  }

  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_velocities[index] *= velocity_damping;
    state._water_heights[index] += state._time_delta * state._water_velocities[index];
  }

  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    float sphere_body_height = 0;
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        sphere_body_height += sphere_body_heights[sphere * state._n * state._m + i * state._m + j];
      }
    }

    const float force = -std::max(sphere_body_height, 0.0f) * std::pow(state._spacing, 2) * GRAVITY;
    const float acceleration = force / sphere_masses[sphere];
    state._sphere_velocities[3 * sphere + 1] += state._time_delta * acceleration;
    // This is what he did in the video, but it doesn't make sense to me.
    if (sphere_body_height > 0)
      state._sphere_velocities[3 * sphere + 1] *= 0.999;
  }

  return state;
}

State apply_sphere_sphere_interaction(State state, const std::vector<float> &sphere_masses, float restitution) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  for (size_t i = 0; i < n_spheres; ++i) {
    for (size_t j = i + 1; j < n_spheres; ++j) {
      const std::array<float, 3> direction = {state._sphere_centers[3 * i] - state._sphere_centers[3 * j],
                                              state._sphere_centers[3 * i + 1] - state._sphere_centers[3 * j + 1],
                                              state._sphere_centers[3 * i + 2] - state._sphere_centers[3 * j + 2]};
      float distance = std::hypot(direction[0], direction[1], direction[2]);
      if (distance <= state._sphere_radii[i] + state._sphere_radii[j]) {
        const std::array<float, 3> unit_direction = {direction[0] / distance, direction[1] / distance,
                                                     direction[2] / distance};

        const float correction = (state._sphere_radii[i] + state._sphere_radii[j] - distance) / 2;
        state._sphere_centers[3 * i] += correction * unit_direction[0];
        if (state._sphere_centers[3 * i + 1] < state._sphere_radii[i])
          state._sphere_centers[3 * i + 1] = state._sphere_radii[i];
        else
          state._sphere_centers[3 * i + 1] += correction * unit_direction[1];
        state._sphere_centers[3 * i + 2] += correction * unit_direction[2];

        state._sphere_centers[3 * j] -= correction * unit_direction[0];
        if (state._sphere_centers[3 * j + 1] < state._sphere_radii[j])
          state._sphere_centers[3 * j + 1] = state._sphere_radii[j];
        else
          state._sphere_centers[3 * j + 1] += correction * unit_direction[1];
        state._sphere_centers[3 * j + 2] -= correction * unit_direction[2];

        const float vi = state._sphere_velocities[3 * i] * unit_direction[0] +
                         state._sphere_velocities[3 * i + 1] * unit_direction[1] +
                         state._sphere_velocities[3 * i + 2] * unit_direction[2];
        const float vj = state._sphere_velocities[3 * j] * unit_direction[0] +
                         state._sphere_velocities[3 * j + 1] * unit_direction[1] +
                         state._sphere_velocities[3 * j + 2] * unit_direction[2];
        const float total_mass = sphere_masses[i] + sphere_masses[j];

        const float new_vi =
            (sphere_masses[i] * vi + sphere_masses[j] * vj - sphere_masses[j] * (vi - vj) * restitution) / total_mass;
        const float new_vj =
            (sphere_masses[i] * vi + sphere_masses[j] * vj - sphere_masses[i] * (vi - vj) * restitution) / total_mass;

        state._sphere_velocities[3 * i] += (new_vi - vi) * unit_direction[0];
        state._sphere_velocities[3 * i + 1] += (new_vi - vi) * unit_direction[1];
        state._sphere_velocities[3 * i + 2] += (new_vi - vi) * unit_direction[2];

        state._sphere_velocities[3 * j] += (new_vj - vj) * unit_direction[0];
        state._sphere_velocities[3 * j + 1] += (new_vj - vj) * unit_direction[1];
        state._sphere_velocities[3 * j + 2] += (new_vj - vj) * unit_direction[2];
      }
    }
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

  State new_state = state;

  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    new_state._sphere_velocities[3 * sphere + 1] += new_state._time_delta * GRAVITY;
    new_state._sphere_centers[3 * sphere + 1] += new_state._time_delta * new_state._sphere_velocities[3 * sphere + 1];
  }

  auto smoothed_sphere_body_heights =
      smooth_sphere_body_heights(sphere_body_heights(new_state._sphere_centers, new_state._sphere_radii,
                                                     new_state._water_xzs, new_state._water_heights),
                                 n_spheres, new_state._n, new_state._m);

  new_state = apply_sphere_water_interaction(new_state, smoothed_sphere_body_heights, sphere_masses);

  constexpr float restitution = 0.1;
  new_state = apply_container_collisions(new_state, restitution);
  new_state = apply_sphere_sphere_interaction(new_state, sphere_masses, restitution);

  return new_state;
}

} // namespace water_simulator::engine