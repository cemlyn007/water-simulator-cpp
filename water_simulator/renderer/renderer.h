#pragma once
#include "water_simulator/renderer/entities/light.h"
#include "water_simulator/renderer/shader.h"
#include <GLFW/glfw3.h>

namespace water_simulator::renderer {

void init();

void terminate();

class Renderer {
private:
  GLFWwindow *_window;
  Shader _shader;
  double _scroll_offset;
  double _mouse_position_in_pixels[2];
  double _last_mouse_position_in_pixels[2];
  double _mouse_position_change_in_pixels[2];
  bool _mouse_click;
  bool _escape_pressed;

  float _camera_radians[2];
  int _framebuffer_width, _framebuffer_height;

  std::array<float, 3> _camera_position;

  entities::Light _light;

public:
  Renderer(int window_width, int window_height);
  ~Renderer();

  void render();
  bool should_close();

private:
  void on_aspect_change();
  void update_camera();

  GLFWwindow *create_window(int width, int height);
  static void key_callback(GLFWwindow *window, int key, int scancode,
                           int action, int mods);
  static void mouse_button_callback(GLFWwindow *window, int button, int action,
                                    int mods);
  static void cursor_position_callback(GLFWwindow *window, double xpos,
                                       double ypos);
  static void scroll_callback(GLFWwindow *window, double xoffset,
                              double yoffset);
  static void framebuffer_size_callback(GLFWwindow *window, int width,
                                        int height);
};

} // namespace water_simulator::renderer
