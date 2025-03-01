#include "water_simulator/engine/state.h"

namespace water_simulator::engine {

constexpr double M_PI = 3.14159265358979323846;

float sphere_mass(float radius, float density) { return 4.0 / 3.0 * M_PI * radius * radius * radius * density; }

State::State(size_t spheres, size_t n, size_t m, float spacing, float water_height, std::vector<float> sphere_radii,
             std::vector<float> sphere_densities)
    : _n(n), _m(m), _spacing(spacing), _sphere_centers(spheres * 3, 0.0), _water_heights(n * m, water_height),
      _water_xzs(), _body_heights(n * m, 0.0), _sphere_velocities(spheres * 3, 0.0), _water_velocities(n * m),
      _wave_speed(2.0), _time_delta(0.0), _sphere_radii(sphere_radii), _sphere_densities(sphere_densities),
      _sphere_masses(spheres), _sphere_body_heights(spheres * n * m, 0.0) {
  for (size_t x = 0; x < _n; ++x) {
    for (size_t z = 0; z < _m; ++z) {
      _water_xzs.push_back(x * _spacing - _spacing * (_n - 1) / 2);
      _water_xzs.push_back(z * _spacing - _spacing * (_m - 1) / 2);
    }
  }
  for (size_t i = 0; i < spheres; ++i)
    _sphere_masses[i] = sphere_mass(_sphere_radii[i], _sphere_densities[i]);
}

} // namespace water_simulator::engine
