#include "math.h"
#include "water_simulator/engine/engine.h"
#include "water_simulator/renderer/renderer.h"
#include <chrono>
#include <iostream>
#include <thread>

using namespace std::chrono_literals;

using namespace water_simulator;

constexpr size_t RESOLUTION = 100;
constexpr float SPACING = 0.02;

int main(int argc, char *argv[]) {
  renderer::init();
  engine::State state(2, RESOLUTION, RESOLUTION, 0.02, 0.8, {0.3, 0.3}, {0.7, 0.7});
  state._sphere_centers[1] = 2.0;
  state._sphere_centers[4] = 5.0;
  renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
                              {{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, 0.3}, {{0.0, 1.0, 0.0}, {0.0, 4.0, 0.0}, 0.3}});
  auto ms = 16ms;
  auto start = std::chrono::high_resolution_clock::now();
  while (!renderer.should_close()) {
    state._time_delta = ms.count() / 1000.0;
    state = engine::step(state);
    renderer.render(state);
    auto end = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Frame time: " << ms.count() << "ms" << std::endl;
    start = end;
    std::this_thread::sleep_for(16ms - ms);
  }
  renderer::terminate();
  return 0;
}