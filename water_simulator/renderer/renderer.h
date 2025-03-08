#pragma once
#include "water_simulator/engine/state.h"
#include "water_simulator/renderer/camera.h"
#include "water_simulator/renderer/entities/ball.h"
#include "water_simulator/renderer/entities/container.h"
#include "water_simulator/renderer/entities/light.h"
#include "water_simulator/renderer/entities/water.h"
#include <GLFW/glfw3.h>

namespace water_simulator::renderer {

struct BallConfig {
  std::array<float, 3> color;
  float radius;
};

void init();

void terminate();

class Renderer {
private:
  GLFWwindow *_window;
  double _scroll_offset;
  double _mouse_position_in_pixels[2];
  double _last_mouse_position_in_pixels[2];
  double _mouse_position_change_in_pixels[2];
  bool _mouse_click;
  bool _escape_pressed;

  float _camera_radians[2];
  int _window_width, _window_height;
  int _framebuffer_width, _framebuffer_height;

  Camera _camera;
  entities::Light _light;
  entities::Container _container;
  entities::Water _water;
  std::vector<entities::Ball> _balls;
  std::array<float, 16> _projection;
  std::array<float, 16> _view;

public:
  Renderer(int window_width, int window_height, size_t resolution, float spacing, float wall_thickness,
           const std::vector<BallConfig> &ball_configs);
  ~Renderer();

  std::array<float, 3> _camera_position;

  void render(const engine::State &state);
  bool should_close();
  std::array<float, 3> get_cursor_direction();

private:
  void on_framebuffer_shape_change();
  void update_camera();

  GLFWwindow *create_window(int width, int height);
  static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods);
  static void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);
  static void cursor_position_callback(GLFWwindow *window, double xpos, double ypos);
  static void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);
  static void framebuffer_size_callback(GLFWwindow *window, int width, int height);
};

} // namespace water_simulator::renderer
