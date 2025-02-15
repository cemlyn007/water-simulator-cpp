#pragma once
#include "water_simulator/engine/state.h"
#include <span>
#include <vector>

namespace water_simulator::engine {

float sphere_mass(float radius, float density);

const std::vector<float> sphere_body_heights(
    const std::vector<float> &sphere_centers,
    const std::vector<float> &sphere_radii, const std::vector<float> &water_xzs,
    const std::vector<float> &water_heights, size_t input_padding);

const std::vector<float> cross_correlation(const std::span<float> &input,
                                           const std::vector<float> &kernel,
                                           const size_t input_n,
                                           const size_t input_m,
                                           const size_t kernel_n,
                                           const size_t kernel_m);

State apply_container_collisions(State state);

State step(const State &state);

} // namespace water_simulator::engine