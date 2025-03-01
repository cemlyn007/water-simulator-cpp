#include "math.h"
#include "water_simulator/engine/engine.h"
#include "water_simulator/renderer/renderer.h"
#include <chrono>
#include <iostream>

using namespace std::chrono_literals;

using namespace water_simulator;

constexpr size_t RESOLUTION = 101;
constexpr float SPACING = 0.02;

int main(int argc, char *argv[]) {
  renderer::init();

  engine::State state(3, RESOLUTION, RESOLUTION, SPACING, 0.8, {0.25, 0.25, 0.25}, {0.265, 0.265, 0.265});

  renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
                              {{{1.0, 0.0, 0.0}, 0.25}, {{0.0, 1.0, 0.0}, 0.25}, {{0.0, 0.0, 1.0}, 0.25}});

  state._sphere_centers[0] = -0.5;
  state._sphere_centers[1] = 3.0;
  state._sphere_centers[2] = -0.5;

  state._sphere_centers[3] = 0.00;
  state._sphere_centers[4] = 3.0;
  state._sphere_centers[5] = 0.00;

  state._sphere_centers[6] = 0.5;
  state._sphere_centers[7] = 3.0;
  state._sphere_centers[8] = 0.5;

  auto us = 1us;
  auto start = std::chrono::high_resolution_clock::now();
  while (!renderer.should_close()) {
    state._time_delta = us.count() / 1000000.0;
    engine::step(state);
    renderer.render(state);
    auto end = std::chrono::high_resolution_clock::now();
    us = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "Frame time: " << us.count() << "us\n";
    std::cout << "Sphere Heights: " << state._sphere_centers[1] << ", " << state._sphere_centers[4] << ", "
              << state._sphere_centers[7] << std::endl;
    start = end;
  }
  renderer::terminate();
  return 0;
}