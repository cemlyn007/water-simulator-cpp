#include "water_simulator/engine/state.h"

namespace water_simulator::engine {
State::State(size_t spheres, size_t n, size_t m, float spacing,
             float water_height, std::vector<float> sphere_radii,
             std::vector<float> sphere_densities)
    : _n(n), _m(m), _spacing(spacing), _sphere_centers(spheres * 3, 0.0),
      _water_heights(n * m, water_height), _water_xzs(),
      _body_heights(n * m, 0.0), _sphere_velocities(spheres * 3, 0.0),
      _water_velocities(n * m), _wave_speed(2.0), _time_delta(0.0),
      _sphere_radii(sphere_radii), _sphere_densities(sphere_densities) {

  // TODO: Initialise water XZs, don't forget to center?
  for (size_t x = 0; x < _n; ++x) {
    for (size_t z = 0; z < _m; ++z) {
      _water_xzs.push_back(x * _spacing);
      _water_xzs.push_back(z * _spacing);
    }
  }
}

State::~State() {}

} // namespace water_simulator::engine
