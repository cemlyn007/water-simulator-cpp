#include "math.h"
#include "water_simulator/engine/engine.h"
#include "water_simulator/renderer/algebra.h"
#include "water_simulator/renderer/renderer.h"
#include <chrono>
#include <iostream>

using namespace std::chrono_literals;

using namespace water_simulator;

constexpr size_t RESOLUTION = 101;
constexpr float SPACING = 0.02;
constexpr float WALL_THICKNESS = 0.01;
constexpr float WALL_HEIGHT = 1.25;

int main(int argc, char *argv[]) {
  renderer::init();

  engine::State state(3, RESOLUTION, RESOLUTION, SPACING, 0.8, {0.25, 0.25, 0.25}, {0.265, 0.265, 0.265});

  renderer::Renderer renderer(1080, 1080, RESOLUTION, SPACING, WALL_THICKNESS,
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
  std::optional<std::pair<size_t, float>> selected_sphere;
  std::array<float, 3> last_position;
  bool rotate_camera = false;
  bool just_selected = false;
  while (!renderer.should_close()) {
    state._time_delta = us.count() / 1000000.0;
    if (rotate_camera) {
      rotate_camera = renderer._mouse_click;
    } else {
      if (renderer._mouse_click) {
        std::array<float, 3> cursor_direction = renderer.get_cursor_direction();
        if (!selected_sphere.has_value()) {
          selected_sphere = engine::raycast(state._sphere_centers, state._sphere_radii, renderer._camera_position,
                                            cursor_direction, SPACING * (RESOLUTION - 1), WALL_THICKNESS, WALL_HEIGHT);
          just_selected = selected_sphere.has_value();
        }
        if (selected_sphere.has_value()) {
          std::cout << "Selected Sphere Index: " << selected_sphere.value().first << std::endl;
          const size_t sphere_index = selected_sphere.value().first;
          const float distance = std::sqrt(selected_sphere.value().second); // perf
          cursor_direction = renderer::normalize(cursor_direction);

          if (just_selected) {
            just_selected = false;
            rotate_camera = false;
            state._sphere_velocities[sphere_index * 3] = 0;
            state._sphere_velocities[sphere_index * 3 + 1] = 0;
            state._sphere_velocities[sphere_index * 3 + 2] = 0;
            last_position[0] = renderer._camera_position[0] + distance * cursor_direction[0];
            last_position[1] = renderer._camera_position[1] + distance * cursor_direction[1];
            last_position[2] = renderer._camera_position[2] + distance * cursor_direction[2];
          } else {
            state._sphere_centers[sphere_index * 3] = renderer._camera_position[0] + distance * cursor_direction[0];
            state._sphere_centers[sphere_index * 3 + 1] = renderer._camera_position[1] + distance * cursor_direction[1];
            state._sphere_centers[sphere_index * 3 + 2] = renderer._camera_position[2] + distance * cursor_direction[2];
            state._sphere_velocities[sphere_index * 3] =
                (state._sphere_centers[sphere_index * 3] - last_position[0]) / state._time_delta;
            state._sphere_velocities[sphere_index * 3 + 1] =
                (state._sphere_centers[sphere_index * 3 + 1] - last_position[1]) / state._time_delta;
            state._sphere_velocities[sphere_index * 3 + 2] =
                (state._sphere_centers[sphere_index * 3 + 2] - last_position[2]) / state._time_delta;
            last_position[0] = state._sphere_centers[sphere_index * 3];
            last_position[1] = state._sphere_centers[sphere_index * 3 + 1];
            last_position[2] = state._sphere_centers[sphere_index * 3 + 2];
          }
        } else {
          rotate_camera = true;
        }
      } else {
        rotate_camera = false;
        selected_sphere.reset();
      }
    }
    engine::step(state);
    renderer.render(state, rotate_camera);
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