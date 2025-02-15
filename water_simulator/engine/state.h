#pragma once
#include <vector>

namespace water_simulator::engine {
class State {
public:
  State(size_t spheres, size_t n, size_t m, float spacing, float water_height,
        std::vector<float> sphere_radii, std::vector<float> sphere_densities);
  ~State();

  size_t _n;
  size_t _m;
  float _spacing;

  std::vector<float> _sphere_centers;
  std::vector<float> _water_heights;
  std::vector<float> _water_xzs;
  std::vector<float> _body_heights;

  std::vector<float> _sphere_velocities;
  std::vector<float> _water_velocities;

  float _wave_speed;
  float _time_delta;

  std::vector<float> _sphere_radii;
  std::vector<float> _sphere_densities;
};
} // namespace water_simulator::engine