#include "math.h"
#include "water_simulator/engine/engine.h"
#include "water_simulator/renderer/renderer.h"
#include <chrono>
#include <iostream>

using namespace std::chrono_literals;

using namespace water_simulator;

constexpr size_t RESOLUTION = 100;
constexpr float SPACING = 0.02;

int main(int argc, char *argv[]) {
  renderer::init();
  engine::State state(3, RESOLUTION, RESOLUTION, SPACING, 0.8, {0.2, 0.3, 0.25}, {2.0, 0.7, 0.1});

  renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
                              {{{1.0, 0.0, 0.0}, 0.2}, {{0.0, 1.0, 0.0}, 0.3}, {{0.0, 0.0, 1.0}, 0.25}});

  state._sphere_centers[0] = -0.5;
  state._sphere_centers[1] = 1.0;
  state._sphere_centers[2] = -0.5;

  state._sphere_centers[3] = 0.00;
  state._sphere_centers[4] = 1.0;
  state._sphere_centers[5] = 0.00;

  state._sphere_centers[6] = 0.5;
  state._sphere_centers[7] = 1.0;
  state._sphere_centers[8] = 0.5;

  auto ms = 1ms;
  auto start = std::chrono::high_resolution_clock::now();
  while (!renderer.should_close()) {
    state._time_delta = ms.count() / 1000.0;
    state = engine::step(state);
    renderer.render(state);
    auto end = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Frame time: " << ms.count() << "ms\n";
    start = end;
  }
  renderer::terminate();
  return 0;
}