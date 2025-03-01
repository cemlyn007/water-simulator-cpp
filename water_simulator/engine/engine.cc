#include "water_simulator/engine/engine.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <mdspan>
#include <vector>

namespace water_simulator::engine {

void sphere_body_heights(std::vector<float> &body_heights, const std::vector<float> &sphere_centers,
                         const std::vector<float> &sphere_radii, const std::vector<float> &water_xzs,
                         const std::vector<float> &water_heights) {
  assert(sphere_centers.size() % 3 == 0);
  assert(water_xzs.size() % 2 == 0);
  assert(sphere_radii.size() * 3 == sphere_centers.size());

  const size_t n_spheres = sphere_centers.size() / 3;
  const size_t n_water_points = water_xzs.size() / 2;

  std::fill(body_heights.begin(), body_heights.end(), 0.0);
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    const double sphere_y = sphere_centers[3 * sphere + 1];
    const double sphere_x = sphere_centers[3 * sphere];
    const double sphere_z = sphere_centers[3 * sphere + 2];
    const double sphere_radius = sphere_radii[sphere];
    const double sphere_radius2 = sphere_radius * sphere_radius;

    for (size_t water_index = 0; water_index < n_water_points; ++water_index) {
      const double x = water_xzs[2 * water_index];
      const double z = water_xzs[2 * water_index + 1];
      const double y = water_heights[water_index];
      const double dx = x - sphere_x;
      const double dz = z - sphere_z;
      const double distance2 = dx * dx + dz * dz;
      if (distance2 < sphere_radius2) {
        const double half_body_height = std::sqrt(sphere_radius2 - distance2);
        const double min_body_height = std::max(sphere_y - half_body_height, 0.0);
        const double max_body_height = std::min(sphere_y + half_body_height, y);
        const double body_height = max_body_height - min_body_height;
        if (body_height > 0) {
          body_heights[sphere * n_water_points + water_index] = body_height;
        }
      }
    }
  }
}

void cross_correlation(std::span<float> output, const std::span<float> input, const std::array<float, 9> &kernel,
                       const size_t input_n, const size_t input_m) {
  constexpr size_t kernel_n = 3;
  constexpr size_t kernel_m = 3;
  if (input.size() != input_n * input_m)
    throw std::runtime_error("Invalid input size");
  if (output.size() != input.size())
    throw std::runtime_error("Invalid output size");
  for (size_t i = 0; i < input_n; ++i) {
    for (size_t j = 0; j < input_m; ++j) {
      float output_element = 0.0;
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
          output_element += input[get_i * input_m + get_j] * kernel[ki * kernel_m + kj];
        }
      }
      output[i * input_m + j] = output_element;
    }
  }
}

void cross_correlation(std::vector<float> &output, const std::span<float> input, const std::array<float, 9> &kernel,
                       const size_t input_n, const size_t input_m) {
  constexpr size_t kernel_n = 3;
  constexpr size_t kernel_m = 3;
  if (input.size() != input_n * input_m)
    throw std::runtime_error("Invalid input size");
  if (output.size() != input.size())
    throw std::runtime_error("Invalid output size");
  for (size_t i = 0; i < input_n; ++i) {
    for (size_t j = 0; j < input_m; ++j) {
      float output_element = 0.0;
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
          output_element += input[get_i * input_m + get_j] * kernel[ki * kernel_m + kj];
        }
      }
      output[i * input_m + j] = output_element;
    }
  }
}

static constexpr std::array<float, 9> SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL{1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f,
                                                                        1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f,
                                                                        1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f};
void smooth_body_heights(std::span<float> body_heights, const size_t n_spheres, const size_t n, const size_t m) {
  std::vector<float> smoothed(n * m);
  for (size_t sphere_index = 0; sphere_index < n_spheres; ++sphere_index) {
    auto sphere_span = body_heights.subspan(sphere_index * n * m, n * m);
    cross_correlation(smoothed, sphere_span, SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL, n, m);
    cross_correlation(sphere_span, std::span(smoothed), SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL, n, m);
  }
}

void apply_container_collisions(State &state, float restitution) {
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
};

static constexpr float GRAVITY = -9.81;
static constexpr std::array<float, 9> NEIGHBOUR_KERNEL{
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
};
static constexpr float NEIGHBOUR_KERNEL_SUM = 4;
void apply_sphere_water_interaction(State &state) {
  const size_t n_spheres = state._sphere_centers.size() / 3;

  assert(state._water_heights.size() == state._n * state._m);
  assert(state._sphere_body_heights.size() == n_spheres * state._n * state._m);
  assert(state._water_velocities.size() == state._water_heights.size());
  assert(state._sphere_velocities.size() == 3 * n_spheres);

  std::vector<float> body_heights(state._n * state._m, 0.0);
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        body_heights[i * state._m + j] += state._sphere_body_heights[sphere * state._n * state._m + i * state._m + j];
      }
    }
  }

  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_heights[index] += 0.05 * (body_heights[index] - state._body_heights[index]);
  }
  state._body_heights = body_heights;

  std::vector<float> sums(state._water_heights.size());
  cross_correlation(sums, state._water_heights, NEIGHBOUR_KERNEL, state._n, state._m);

  double wave_speed = std::min(state._wave_speed, 0.5 * state._spacing / state._time_delta);
  const float c = std::pow(wave_speed / state._spacing, 2);

  constexpr double POSITIONAL_DAMPING = 1.0f;
  double positional_damping = std::min(POSITIONAL_DAMPING * state._time_delta, 1.0);
  constexpr float VELOCITY_DAMPING = 0.3;
  auto velocity_damping = std::max(0.0, 1.0 - VELOCITY_DAMPING * state._time_delta);
  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_velocities[index] +=
        state._time_delta * c * (sums[index] - NEIGHBOUR_KERNEL_SUM * state._water_heights[index]);
    state._water_heights[index] +=
        (sums[index] / NEIGHBOUR_KERNEL_SUM - state._water_heights[index]) * positional_damping;
  }

  for (size_t index = 0; index < state._water_heights.size(); ++index) {
    state._water_velocities[index] *= velocity_damping;
    state._water_heights[index] += state._time_delta * state._water_velocities[index];
  }

  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    float sphere_body_height = 0;
    for (size_t i = 0; i < state._n; ++i) {
      for (size_t j = 0; j < state._m; ++j) {
        sphere_body_height += state._sphere_body_heights[sphere * state._n * state._m + i * state._m + j];
      }
    }

    const float force = -std::max(sphere_body_height, 0.0f) * std::pow(state._spacing, 2) * GRAVITY;
    const float acceleration = force / state._sphere_masses[sphere];

    if (sphere_body_height > 0) {
      state._sphere_velocities[3 * sphere + 1] += state._time_delta * acceleration;
      state._sphere_velocities[3 * sphere + 1] *= 0.999;
    }
  }
}

void apply_sphere_sphere_interaction(std::vector<float> &centers, std::vector<float> &velocities,
                                     const std::vector<float> &radii, const std::vector<float> &masses,
                                     float restitution) {
  const size_t n_spheres = centers.size() / 3;
  for (size_t i = 0; i < n_spheres; ++i) {
    for (size_t j = i + 1; j < n_spheres; ++j) {
      const std::array<float, 3> direction = {centers[3 * i] - centers[3 * j], centers[3 * i + 1] - centers[3 * j + 1],
                                              centers[3 * i + 2] - centers[3 * j + 2]};
      float distance = std::hypot(direction[0], direction[1], direction[2]);
      if (distance <= radii[i] + radii[j]) {
        const std::array<float, 3> unit_direction = {direction[0] / distance, direction[1] / distance,
                                                     direction[2] / distance};

        const float correction = (radii[i] + radii[j] - distance) / 2;
        centers[3 * i] += correction * unit_direction[0];
        if (centers[3 * i + 1] < radii[i])
          centers[3 * i + 1] = radii[i];
        else
          centers[3 * i + 1] += correction * unit_direction[1];
        centers[3 * i + 2] += correction * unit_direction[2];

        centers[3 * j] -= correction * unit_direction[0];
        if (centers[3 * j + 1] < radii[j])
          centers[3 * j + 1] = radii[j];
        else
          centers[3 * j + 1] += correction * unit_direction[1];
        centers[3 * j + 2] -= correction * unit_direction[2];

        const float vi = velocities[3 * i] * unit_direction[0] + velocities[3 * i + 1] * unit_direction[1] +
                         velocities[3 * i + 2] * unit_direction[2];
        const float vj = velocities[3 * j] * unit_direction[0] + velocities[3 * j + 1] * unit_direction[1] +
                         velocities[3 * j + 2] * unit_direction[2];
        const float total_mass = masses[i] + masses[j];

        const float new_vi = (masses[i] * vi + masses[j] * vj - masses[j] * (vi - vj) * restitution) / total_mass;
        const float new_vj = (masses[i] * vi + masses[j] * vj - masses[i] * (vi - vj) * restitution) / total_mass;

        velocities[3 * i] += (new_vi - vi) * unit_direction[0];
        velocities[3 * i + 1] += (new_vi - vi) * unit_direction[1];
        velocities[3 * i + 2] += (new_vi - vi) * unit_direction[2];

        velocities[3 * j] += (new_vj - vj) * unit_direction[0];
        velocities[3 * j + 1] += (new_vj - vj) * unit_direction[1];
        velocities[3 * j + 2] += (new_vj - vj) * unit_direction[2];
      }
    }
  }
}

void step(State &state) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  if (n_spheres != state._sphere_velocities.size() / 3 || n_spheres != state._sphere_radii.size() ||
      n_spheres != state._sphere_densities.size()) {
    throw std::runtime_error("Invalid state: sphere vectors have different sizes");
  }

  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    state._sphere_velocities[3 * sphere + 1] += state._time_delta * GRAVITY;
    state._sphere_centers[3 * sphere + 1] += state._time_delta * state._sphere_velocities[3 * sphere + 1];
  }

  sphere_body_heights(state._sphere_body_heights, state._sphere_centers, state._sphere_radii, state._water_xzs,
                      state._water_heights);
  smooth_body_heights(state._sphere_body_heights, n_spheres, state._n, state._m);

  apply_sphere_water_interaction(state);

  constexpr float RESTITUTION = 0.1;
  apply_container_collisions(state, RESTITUTION);
  apply_sphere_sphere_interaction(state._sphere_centers, state._sphere_velocities, state._sphere_radii,
                                  state._sphere_masses, RESTITUTION);
}

} // namespace water_simulator::engine