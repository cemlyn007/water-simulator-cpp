#include "water_simulator/engine/engine.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <sycl/sycl.hpp>
#include <vector>

namespace water_simulator::engine {

static sycl::queue queue;

sycl::event sphere_body_heights(sycl::buffer<float, 1> &d_body_heights, sycl::buffer<float, 1> &d_sphere_centers,
                                sycl::buffer<float, 1> &d_sphere_radii, sycl::buffer<float, 1> &d_water_xzs,
                                sycl::buffer<float, 1> &d_water_heights) {
  assert(sphere_centers.size() % 3 == 0);
  assert(water_xzs.size() % 2 == 0);
  assert(sphere_radii.size() * 3 == sphere_centers.size());

  const size_t n_spheres = d_sphere_centers.size() / 3;
  const size_t n_water_points = d_water_xzs.size() / 2;
  return queue.submit([&](sycl::handler &handler) {
    auto sphere_centers_acc = d_sphere_centers.get_access<sycl::access::mode::read>(handler);
    auto sphere_radii_acc = d_sphere_radii.get_access<sycl::access::mode::read>(handler);
    auto water_xzs_acc = d_water_xzs.get_access<sycl::access::mode::read>(handler);
    auto water_heights_acc = d_water_heights.get_access<sycl::access::mode::read>(handler);
    auto body_heights_acc = d_body_heights.get_access<sycl::access::mode::discard_write>(handler);

    handler.parallel_for(sycl::range<1>(n_water_points), [=](sycl::id<1> index) {
      auto water_index = index;

      const float x = water_xzs_acc[2 * water_index];
      const float z = water_xzs_acc[2 * water_index + 1];
      const float y = water_heights_acc[water_index];

      for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
        const float sphere_y = sphere_centers_acc[3 * sphere + 1];
        const float sphere_x = sphere_centers_acc[3 * sphere];
        const float sphere_z = sphere_centers_acc[3 * sphere + 2];
        const float sphere_radius = sphere_radii_acc[sphere];
        const float sphere_radius2 = sphere_radius * sphere_radius;
        const float dx = x - sphere_x;
        const float dz = z - sphere_z;
        const float distance2 = dx * dx + dz * dz;

        body_heights_acc[sphere * n_water_points + water_index] = 0;
        if (distance2 < sphere_radius2) {
          const float half_body_height = std::sqrt(sphere_radius2 - distance2);
          const float min_body_height = std::max(sphere_y - half_body_height, 0.0f);
          const float max_body_height = std::min(sphere_y + half_body_height, y);
          const float body_height = max_body_height - min_body_height;
          body_heights_acc[sphere * n_water_points + water_index] = std::max(body_height, 0.0f);
        }
      }
    });
  });
}

template <typename S, typename T>
sycl::event cross_correlation_impl(sycl::event dependency, S &d_output, T &d_input, sycl::buffer<float, 1> &d_kernel,
                                   const size_t input_n, const size_t input_m) {
  constexpr size_t kernel_n = 3;
  constexpr size_t kernel_m = 3;
  if (d_input.size() != input_n * input_m)
    throw std::runtime_error("Invalid input size");
  if (d_output.size() != d_input.size())
    throw std::runtime_error("Invalid output size");
  return queue.submit([&](sycl::handler &handler) {
    handler.depends_on(dependency);
    auto input_acc = d_input.template get_access<sycl::access::mode::read>(handler);
    auto kernel_acc = d_kernel.template get_access<sycl::access::mode::read>(handler);
    auto output_acc = d_output.template get_access<sycl::access::mode::discard_write>(handler);

    handler.parallel_for(sycl::range<2>({input_n, input_m}), [=](sycl::id<2> index) {
      size_t i = index[0];
      size_t j = index[1];

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
          output_element += input_acc[get_i * input_m + get_j] * kernel_acc[ki * kernel_m + kj];
        }
      }
      output_acc[i * input_m + j] = output_element;
    });
  });
  ;
}

template <typename S, typename T>
sycl::event cross_correlation_impl(sycl::event dependency, S &d_output, T &d_input, sycl::buffer<float, 1> &d_kernel,
                                   const size_t input_n, const size_t input_m, const size_t input_k) {
  constexpr size_t kernel_n = 3;
  constexpr size_t kernel_m = 3;
  if (d_input.size() != input_n * input_m * input_k)
    throw std::runtime_error("Invalid input size");
  if (d_output.size() != d_input.size())
    throw std::runtime_error("Invalid output size");
  return queue.submit([&](sycl::handler &handler) {
    handler.depends_on(dependency);
    auto input_acc = d_input.template get_access<sycl::access::mode::read>(handler);
    auto kernel_acc = d_kernel.template get_access<sycl::access::mode::read>(handler);
    auto output_acc = d_output.template get_access<sycl::access::mode::discard_write>(handler);
    handler.parallel_for(sycl::range<3>({input_n, input_m, input_k}), [=](sycl::id<3> index) {
      size_t i = index[0];
      size_t j = index[1];
      size_t k = index[2];
      size_t offset = i * input_m * input_k;
      float output_element = 0.0;
      for (size_t kj = 0; kj < kernel_n; ++kj) {
        size_t get_j = j + kj - kernel_n / 2;
        if (j + kj < kernel_n / 2) {
          get_j = 0;
        } else if (get_j >= input_m) {
          get_j = input_m - 1;
        }
        for (size_t kk = 0; kk < kernel_m; ++kk) {
          size_t get_k = k + kk - kernel_m / 2;
          if (k + kk < kernel_m / 2) {
            get_k = 0;
          } else if (get_k >= input_k) {
            get_k = input_k - 1;
          }
          output_element += input_acc[offset + get_j * input_k + get_k] * kernel_acc[kj * kernel_m + kk];
        }
      }
      output_acc[offset + j * input_k + k] = output_element;
    });
  });
}

sycl::event cross_correlation(sycl::event dependency, sycl::buffer<float, 1> &d_output, sycl::buffer<float, 1> &d_input,
                              sycl::buffer<float, 1> &d_kernel, const size_t input_n, const size_t input_m) {
  return cross_correlation_impl<sycl::buffer<float, 1> &, sycl::buffer<float, 1> &>(dependency, d_output, d_input,
                                                                                    d_kernel, input_n, input_m);
}

sycl::event cross_correlation(sycl::event dependency, sycl::buffer<float, 1> &d_output, sycl::buffer<float, 1> &d_input,
                              sycl::buffer<float, 1> &d_kernel, const size_t input_n, const size_t input_m,
                              const size_t input_k) {
  return cross_correlation_impl<sycl::buffer<float, 1> &, sycl::buffer<float, 1> &>(
      dependency, d_output, d_input, d_kernel, input_n, input_m, input_k);
}

static constexpr std::array<float, 9> SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL{1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f,
                                                                        1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f,
                                                                        1.0f / 9.0f, 1.0f / 9.0f, 1.0f / 9.0f};

sycl::event smooth_body_heights(sycl::event dependency, sycl::buffer<float, 1> &d_body_heights,
                                sycl::buffer<float, 1> &d_smoothed, const size_t n_spheres, const size_t n,
                                const size_t m) {
  auto d_kernel = sycl::buffer<float, 1>(SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL.data(),
                                         sycl::range<1>(SMOOTH_SPHERE_BODY_HEIGHTS_KERNEL.size()));
  auto first_dependency = cross_correlation(dependency, d_smoothed, d_body_heights, d_kernel, n_spheres, n, m);
  return cross_correlation(first_dependency, d_body_heights, d_smoothed, d_kernel, n_spheres, n, m);
}

void apply_container_collisions(State &state, float restitution) {
  const size_t n_spheres = state._sphere_centers.size() / 3;
  for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
    size_t sphere_y_index = 3 * sphere + 1;
    if (state._sphere_centers[sphere_y_index] <= state._sphere_radii[sphere]) {
      state._sphere_centers[sphere_y_index] = state._sphere_radii[sphere];
      state._sphere_velocities[sphere_y_index] = -state._sphere_velocities[sphere_y_index] * restitution;
    }

    const float half_wall_size_n = ((state._n - 1) * state._spacing) / 2.0f;
    size_t sphere_z_index = sphere * 3 + 2;
    // North wall:
    if (half_wall_size_n - state._sphere_centers[sphere_z_index] < state._sphere_radii[sphere]) {
      state._sphere_centers[sphere_z_index] = half_wall_size_n - state._sphere_radii[sphere];
      state._sphere_velocities[sphere_z_index] = -state._sphere_velocities[sphere_z_index] * restitution;
    }
    // South wall:
    if (state._sphere_centers[sphere_z_index] - state._sphere_radii[sphere] < -half_wall_size_n) {
      state._sphere_centers[sphere_z_index] = -half_wall_size_n + state._sphere_radii[sphere];
      state._sphere_velocities[sphere_z_index] = -state._sphere_velocities[sphere_z_index] * restitution;
    }

    const float half_wall_size_m = ((state._m - 1) * state._spacing) / 2.0f;
    size_t sphere_x_index = sphere * 3;
    // East wall:
    if (half_wall_size_m - state._sphere_centers[sphere_x_index] < state._sphere_radii[sphere]) {
      state._sphere_centers[sphere_x_index] = half_wall_size_m - state._sphere_radii[sphere];
      state._sphere_velocities[sphere_x_index] = -state._sphere_velocities[sphere_x_index] * restitution;
    }
    // West wall:
    if (state._sphere_centers[sphere_x_index] - state._sphere_radii[sphere] < -half_wall_size_m) {
      state._sphere_centers[sphere_x_index] = -half_wall_size_m + state._sphere_radii[sphere];
      state._sphere_velocities[sphere_x_index] = -state._sphere_velocities[sphere_x_index] * restitution;
    }
  }
};

static constexpr float GRAVITY = -9.81;
static constexpr float ALPHA = 0.5;
sycl::event apply_body_height_change(sycl::event dependency, sycl::buffer<float, 1> &d_water_heights,
                                     sycl::buffer<float, 1> &d_sphere_body_heights,
                                     sycl::buffer<float, 1> &d_sphere_masses,
                                     sycl::buffer<float, 1> &d_sphere_velocities,
                                     sycl::buffer<float, 1> &d_body_heights, const size_t n, const size_t m,
                                     const double spacing, const double time_delta) {
  const size_t n_spheres = d_sphere_masses.size();
  const size_t n_water_points = n * m;
  const float spacing_squared = spacing * spacing;
  return queue.submit([&](sycl::handler &handler) {
    handler.depends_on(dependency);
    auto water_heights_acc = d_water_heights.get_access<sycl::access::mode::read_write>(handler);
    auto sphere_body_heights_acc = d_sphere_body_heights.get_access<sycl::access::mode::read>(handler);
    auto sphere_masses_acc = d_sphere_masses.get_access<sycl::access::mode::read>(handler);
    auto sphere_velocities_acc = d_sphere_velocities.get_access<sycl::access::mode::read_write>(handler);
    auto body_heights_acc = d_body_heights.get_access<sycl::access::mode::read_write>(handler);
    handler.parallel_for(sycl::range<1>(n_water_points), [=](sycl::id<1> index) {
      float body_height = 0;
      for (size_t sphere = 0; sphere < n_spheres; ++sphere) {
        const float sphere_body_height = sphere_body_heights_acc[sphere * n * m + index];
        const float force = -std::max(sphere_body_height, 0.0f) * spacing_squared * GRAVITY;
        const float acceleration = force / sphere_masses_acc[sphere];
        if (sphere_body_height > 0) {
          sphere_velocities_acc[3 * sphere + 1] += time_delta * acceleration;
          sphere_velocities_acc[3 * sphere + 1] *= 0.999;
        }
        body_height += sphere_body_height;
      }
      water_heights_acc[index] += ALPHA * (body_height - body_heights_acc[index]);
      body_heights_acc[index] = body_height;
    });
  });
}

static constexpr std::array<float, 9> NEIGHBOUR_KERNEL{
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
};
static constexpr float NEIGHBOUR_KERNEL_SUM = 4;
sycl::event apply_neighbour_deltas(sycl::event dependency, sycl::buffer<float, 1> &d_water_heights,
                                   sycl::buffer<float, 1> &d_water_velocities, sycl::buffer<float, 1> &d_neighbour_sums,
                                   sycl::buffer<float, 1> &d_neighbour_kernel, const size_t n, const size_t m,
                                   const double time_delta, const double spacing, const double wave_speed) {
  auto cross_correlation_dependency =
      cross_correlation(dependency, d_neighbour_sums, d_water_heights, d_neighbour_kernel, n, m);
  const size_t n_water_points = n * m;
  const double adjusted_wave_speed = std::min(wave_speed, 0.5 * spacing / time_delta);
  const float c = std::pow(adjusted_wave_speed / spacing, 2);
  constexpr double POSITIONAL_DAMPING = 1.0f;
  double positional_damping = std::min(POSITIONAL_DAMPING * time_delta, 1.0);
  return queue.submit([&](sycl::handler &handler) {
    handler.depends_on(cross_correlation_dependency);
    auto water_heights_acc = d_water_heights.get_access<sycl::access::mode::read_write>(handler);
    auto water_velocities_acc = d_water_velocities.get_access<sycl::access::mode::read_write>(handler);
    auto neighbour_sums_acc = d_neighbour_sums.get_access<sycl::access::mode::read>(handler);
    handler.parallel_for(sycl::range<1>(n_water_points), [=](sycl::id<1> index) {
      water_velocities_acc[index] +=
          time_delta * c * (neighbour_sums_acc[index] - NEIGHBOUR_KERNEL_SUM * water_heights_acc[index]);
      water_heights_acc[index] +=
          (neighbour_sums_acc[index] / NEIGHBOUR_KERNEL_SUM - water_heights_acc[index]) * positional_damping;
    });
  });
}

sycl::event apply_velocity_damping(sycl::event dependency, sycl::buffer<float, 1> d_water_heights,
                                   sycl::buffer<float, 1> d_water_velocities, const size_t n, const size_t m,
                                   const double time_delta) {
  constexpr float VELOCITY_DAMPING = 0.3;
  const float velocity_damping = std::max(0.0, 1.0 - VELOCITY_DAMPING * time_delta);
  const size_t n_water_points = n * m;
  return queue.submit([&](sycl::handler &handler) {
    handler.depends_on(dependency);
    auto water_heights_acc = d_water_heights.get_access<sycl::access::mode::read_write>(handler);
    auto water_velocities_acc = d_water_velocities.get_access<sycl::access::mode::read_write>(handler);
    handler.parallel_for(sycl::range<1>(n_water_points), [=](sycl::id<1> index) {
      water_velocities_acc[index] *= velocity_damping;
      water_heights_acc[index] += time_delta * water_velocities_acc[index];
    });
  });
}

sycl::event
apply_sphere_water_interaction(sycl::event dependency, sycl::buffer<float, 1> &d_water_heights,
                               sycl::buffer<float, 1> &d_sphere_body_heights, sycl::buffer<float, 1> &d_sphere_masses,
                               sycl::buffer<float, 1> &d_sphere_velocities, sycl::buffer<float, 1> &d_body_heights,
                               sycl::buffer<float, 1> &d_water_velocities, sycl::buffer<float, 1> &d_neighbour_sums,
                               sycl::buffer<float, 1> &d_neighbour_kernel, const size_t n, const size_t m,
                               const double spacing, const double time_delta, const double wave_speed) {
  assert(state._water_heights.size() == state._n * state._m);
  assert(state._sphere_body_heights.size() == d_sphere_masses.size() * state._n * state._m);
  assert(state._water_velocities.size() == state._water_heights.size());
  assert(state._sphere_velocities.size() == 3 * d_sphere_masses.size());
  auto apply_body_height_change_event =
      apply_body_height_change(dependency, d_water_heights, d_sphere_body_heights, d_sphere_masses, d_sphere_velocities,
                               d_body_heights, n, m, spacing, time_delta);
  auto apply_neighbour_deltas_event =
      apply_neighbour_deltas(apply_body_height_change_event, d_water_heights, d_water_velocities, d_neighbour_sums,
                             d_neighbour_kernel, n, m, time_delta, spacing, wave_speed);
  return apply_velocity_damping(apply_neighbour_deltas_event, d_water_heights, d_water_velocities, n, m, time_delta);
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
  sycl::buffer<float, 1> d_sphere_body_heights(state._sphere_body_heights.data(),
                                               sycl::range<1>(state._sphere_body_heights.size()));
  sycl::buffer<float, 1> d_sphere_centers(state._sphere_centers.data(), sycl::range<1>(state._sphere_centers.size()));
  sycl::buffer<float, 1> d_sphere_radii(state._sphere_radii.data(), sycl::range<1>(state._sphere_radii.size()));
  sycl::buffer<float, 1> d_water_xzs(state._water_xzs.data(), sycl::range<1>(state._water_xzs.size()));
  sycl::buffer<float, 1> d_water_heights(state._water_heights.data(), sycl::range<1>(state._water_heights.size()));
  static std::vector<float> smooth(n_spheres * state._n * state._m, 0);
  sycl::buffer<float, 1> d_smooth(smooth.data(), sycl::range<1>(smooth.size()));
  sycl::buffer<float, 1> d_sphere_masses(state._sphere_masses.data(), sycl::range<1>(state._sphere_masses.size()));
  sycl::buffer<float, 1> d_sphere_velocities(state._sphere_velocities.data(),
                                             sycl::range<1>(state._sphere_velocities.size()));
  sycl::buffer<float, 1> d_body_heights(state._body_heights.data(), sycl::range<1>(state._body_heights.size()));
  sycl::buffer<float, 1> d_water_velocities(state._water_velocities.data(),
                                            sycl::range<1>(state._water_velocities.size()));
  sycl::buffer<float, 1> d_neighbour_sums(state._neighbour_sums.data(), sycl::range<1>(state._neighbour_sums.size()));
  sycl::buffer<float, 1> d_neighbour_kernel(NEIGHBOUR_KERNEL.data(), sycl::range<1>(NEIGHBOUR_KERNEL.size()));
  auto sphere_body_heights_event =
      sphere_body_heights(d_sphere_body_heights, d_sphere_centers, d_sphere_radii, d_water_xzs, d_water_heights);
  auto smooth_body_heights_event =
      smooth_body_heights(sphere_body_heights_event, d_sphere_body_heights, d_smooth, n_spheres, state._n, state._m);
  apply_sphere_water_interaction(smooth_body_heights_event, d_water_heights, d_sphere_body_heights, d_sphere_masses,
                                 d_sphere_velocities, d_body_heights, d_water_velocities, d_neighbour_sums,
                                 d_neighbour_kernel, state._n, state._m, state._spacing, state._time_delta,
                                 state._wave_speed)
      .wait();
  constexpr float RESTITUTION = 0.1;
  apply_sphere_sphere_interaction(state._sphere_centers, state._sphere_velocities, state._sphere_radii,
                                  state._sphere_masses, RESTITUTION);
  apply_container_collisions(state, RESTITUTION);
}

} // namespace water_simulator::engine