#pragma once
#include <vector>

namespace water_simulator::engine {
class State {
public:
  State(size_t spheres, size_t n, size_t m, float spacing, float water_height, std::vector<float> sphere_radii,
        std::vector<float> sphere_densities);

  size_t _n;
  size_t _m;
  float _spacing;

  std::vector<float> _sphere_centers;
  std::vector<float> _water_heights;
  std::vector<float> _water_xzs;
  std::vector<float> _body_heights;

  std::vector<float> _sphere_velocities;
  std::vector<float> _water_velocities;

  double _wave_speed;
  double _time_delta;

  std::vector<float> _sphere_radii;
  std::vector<float> _sphere_densities;
  std::vector<float> _sphere_masses;

  std::vector<float> _sphere_body_heights;
};
} // namespace water_simulator::engine