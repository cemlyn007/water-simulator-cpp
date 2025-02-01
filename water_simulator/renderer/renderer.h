#ifndef WATER_SIMULATOR_RENDERER_H_
#define WATER_SIMULATOR_RENDERER_H_
#include "water_simulator/renderer/shader.h"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdlib.h>

namespace water_simulator::renderer {

void init();

void terminate();

class Renderer {
private:
  GLFWwindow *_window;
  Shader _shader;
  double _mouse_position_in_pixels[2];
  bool _mouse_click;
  bool _escape_pressed;

public:
  Renderer(int window_width, int window_height);
  ~Renderer();

  void render();
  bool should_close();
  std::tuple<int, bool> get_selected_location();

private:
  GLFWwindow *create_window(int width, int height);
  static void key_callback(GLFWwindow *window, int key, int scancode,
                           int action, int mods);
  static void mouse_button_callback(GLFWwindow *window, int button, int action,
                                    int mods);
  static void cursor_position_callback(GLFWwindow *window, double xpos,
                                       double ypos);
  std::string read_shader(const std::string &file_path);
};

} // namespace water_simulator::renderer

#endif
