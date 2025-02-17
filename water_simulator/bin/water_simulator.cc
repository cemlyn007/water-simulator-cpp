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
  // engine::State state(3, RESOLUTION, RESOLUTION, 0.02, 0.8, {0.2, 0.3, 0.25}, {2.0, 0.7, 0.1});
  // state._sphere_centers[1] = 20.0;
  // state._sphere_centers[4] = 50.0;
  // state._sphere_centers[6] = -0.5;
  // state._sphere_centers[7] = 1.0;
  // state._sphere_centers[8] = -0.5;
  // renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
  //                             {{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, 0.3},
  //                              {{0.0, 1.0, 0.0}, {0.0, 4.0, 0.0}, 0.3},
  //                              {{0.0, 0.0, 1.0}, {0.6, 1.0, 0.6}, 0.1}});

  // // Test entity permutation shows an issue?
  // engine::State state(3, RESOLUTION, RESOLUTION, 0.02, 0.8, {0.25, 0.3, 0.2}, {0.1, 0.7, 2.0});
  // state._sphere_centers[0] = -0.5;
  // state._sphere_centers[1] = 1.0;
  // state._sphere_centers[2] = -0.5;
  // state._sphere_centers[4] = 50.0;
  // state._sphere_centers[7] = 20.0;
  // renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
  //                             {{{0.0, 0.0, 1.0}, {0.6, 1.0, 0.6}, 0.1},
  //                              {{0.0, 1.0, 0.0}, {0.0, 4.0, 0.0}, 0.3},
  //                              {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, 0.3}});

  // Test entity permutation shows an issue?
  engine::State state(3, RESOLUTION, RESOLUTION, SPACING, 0.8, {0.25, 0.25, 0.25}, {0.1, 0.1, 0.1});

  // TEST 1
  // state._sphere_centers[0] = -0.7;
  // state._sphere_centers[1] = 20.0;
  // state._sphere_centers[2] = -0.7;

  // // state._sphere_centers[3] = -0.5;
  // // state._sphere_centers[4] = 20.0;
  // // state._sphere_centers[5] = 0.5;

  // // // Did one extra bounce!
  // // state._sphere_centers[3] = 0.5;
  // // state._sphere_centers[4] = 20.0;
  // // state._sphere_centers[5] = 0.5;

  // //
  // state._sphere_centers[3] = -0.30;
  // state._sphere_centers[4] = 20.0;
  // state._sphere_centers[5] = -0.30;

  // // // Did one extra bounce!
  // state._sphere_centers[6] = 0.5;
  // state._sphere_centers[7] = 20.0;
  // state._sphere_centers[8] = 0.5;

  // // // Did not bounce!
  // // state._sphere_centers[6] = 0.6;
  // // state._sphere_centers[7] = 20.0;
  // // state._sphere_centers[8] = -0.6;

  // // TODO: If the location causes the extra bounce, then bug?

  // TEST 2

  state._sphere_centers[0] = -0.5;
  state._sphere_centers[1] = 20.0;
  state._sphere_centers[2] = -0.5;

  state._sphere_centers[3] = 0.00;
  state._sphere_centers[4] = 20.0;
  state._sphere_centers[5] = 0.00;

  state._sphere_centers[6] = 0.5;
  state._sphere_centers[7] = 20.0;
  state._sphere_centers[8] = 0.5;

  // TODO: If the location causes the extra bounce, then bug?

  renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING,
                              {{{1.0, 0.0, 0.0}, 0.25}, {{0.0, 1.0, 0.0}, 0.25}, {{0.0, 0.0, 1.0}, 0.25}});

  // engine::State state(1, RESOLUTION, RESOLUTION, 0.02, 0.8,
  //                     {
  //                         0.2,
  //                     },
  //                     {
  //                         2.0,
  //                     });
  // state._sphere_centers[1] = 20.0;
  // // state._sphere_centers[4] = 50.0;
  // // state._sphere_centers[6] = -0.5;
  // // state._sphere_centers[7] = 1.0;
  // // state._sphere_centers[8] = -0.5;
  // renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING, {{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, 0.3}});
  auto ms = 16ms;
  auto start = std::chrono::high_resolution_clock::now();
  while (!renderer.should_close()) {
    // state._time_delta = (512ms).count() / 1000.0f; // ms.count() / 1000.0;

    // state._time_delta = (1ms).count() / 1000.0f; // ms.count() / 1000.0;
    state._time_delta = 1.0 / 30.0;
    std::cout << state._time_delta << std::endl;
    state = engine::step(state);
    // state._sphere_centers[4] = 50.0;
    renderer.render(state);
    auto end = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Frame time: " << ms.count() << "ms" << std::endl;
    start = end;
    // std::this_thread::sleep_for(16ms - ms);
  }
  renderer::terminate();
  return 0;
}