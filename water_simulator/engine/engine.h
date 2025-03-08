#pragma once
#include "water_simulator/engine/state.h"
#include <optional>
#include <vector>
namespace water_simulator::engine {

void step(State &state);

std::optional<size_t> raycast(const std::vector<float> &sphere_centers, const std::vector<float> &sphere_radii,
                              const std::array<float, 3> &ray_start, const std::array<float, 3> &ray_direction,
                              const float wall_size, const float wall_thickness, const float wall_height);

} // namespace water_simulator::engine