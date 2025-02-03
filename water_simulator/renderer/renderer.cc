#include "water_simulator/renderer/renderer.h"
#include "water_simulator/renderer/algebra.h"
#include "water_simulator/renderer/gl_error_macro.h"
#include "water_simulator/renderer/shader.h"
#include <GL/glew.h>
#include <iostream>

static const char *SHADER_VERTEX_FILE_PATH =
    "water_simulator/renderer/shaders/shader.vert";
static const char *SHADER_FRAGMENT_FILE_PATH =
    "water_simulator/renderer/shaders/shader.frag";

static void glfwErrorCallback(int error, const char *description) {
  std::cerr << "GLFW Error: " << error << " - " << description << std::endl;
  throw std::runtime_error(std::string("GLFW Error: ") + description);
}

namespace water_simulator::renderer {

void init() {
  glfwSetErrorCallback(glfwErrorCallback);
  if (!glfwInit()) {
    throw std::runtime_error("Could not initialize glfw");
  }
  // else...
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
};

void terminate() {
  glfwSetErrorCallback(nullptr);
  glfwTerminate();
};

Renderer::Renderer(int window_width, int window_height)
    : _window(create_window(window_width, window_height)),
      _shader(read_file(SHADER_VERTEX_FILE_PATH),
              read_file(SHADER_FRAGMENT_FILE_PATH)),
      _mouse_click(false), _escape_pressed(false), _light() {
  on_aspect_change();
  _light.set_color({1.0, 1.0, 1.0});
  _light.set_model(translate(scale(eye4d(), {0.2, 0.2, 0.2}), {1.2, 4.0, 2.0}));
}

Renderer::~Renderer() { glfwDestroyWindow(_window); };

void Renderer::render() {
  glfwMakeContextCurrent(_window);
  GL_CALL(glClear(GL_COLOR_BUFFER_BIT));
  update_camera();
  _shader.use();
  _light.draw();
  GL_CALL(glfwSwapBuffers(_window));
  GL_CALL(glfwPollEvents());
}

void Renderer::on_aspect_change() {
  float aspect = static_cast<float>(_framebuffer_width) /
                 static_cast<float>(_framebuffer_height);
  _light.set_projection(perspective((M_PI * 60) / 180, aspect, 0.01, 100.0));
};

void Renderer::update_camera() {
  _camera_position = {2.5, 3.54, 2.5};
  _light.set_view(look_at(_camera_position, {1.2, 4.0, 2.0}, {0.0, 1.0, 0.0}));
};

bool Renderer::should_close() {
  return glfwWindowShouldClose(_window) || _escape_pressed;
}

GLFWwindow *Renderer::create_window(int width, int height) {
  GLFWwindow *window =
      glfwCreateWindow(width, height, "Water Simulator", nullptr, nullptr);
  if (!window) {
    throw std::runtime_error("Could not create window");
  }
  // else...
  glfwMakeContextCurrent(window);
  GLenum error = glewInit();
  if (GLEW_OK != error) {
    glfwDestroyWindow(window);
    throw std::runtime_error(
        std::string("Error initializing glew: ") +
        reinterpret_cast<const char *>(glewGetErrorString(error)));
  }
  glfwGetFramebufferSize(window, &_framebuffer_width, &_framebuffer_height);
  glfwSetWindowUserPointer(window, this);
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GLFW_TRUE);
  glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);
  glfwSetKeyCallback(window, key_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  return window;
};

void Renderer::key_callback(GLFWwindow *window, int key, int scancode,
                            int action, int mods) {
  Renderer *renderer =
      static_cast<Renderer *>(glfwGetWindowUserPointer(window));
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    (renderer->_escape_pressed) = true;
}

void Renderer::mouse_button_callback(GLFWwindow *window, int button, int action,
                                     int mods) {
  Renderer *renderer =
      static_cast<Renderer *>(glfwGetWindowUserPointer(window));
  if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    (renderer->_mouse_click) = true;
}

void Renderer::cursor_position_callback(GLFWwindow *window, double xpos,
                                        double ypos) {
  Renderer *renderer =
      static_cast<Renderer *>(glfwGetWindowUserPointer(window));
  (renderer->_mouse_position_in_pixels[0]) = xpos;
  (renderer->_mouse_position_in_pixels[1]) = ypos;
}

void Renderer::framebuffer_size_callback(GLFWwindow *window, int width,
                                         int height) {
  Renderer *renderer =
      static_cast<Renderer *>(glfwGetWindowUserPointer(window));

  glfwMakeContextCurrent(window);
  glViewport(0, 0, width, height);
  glfwGetFramebufferSize(window, &(renderer->_framebuffer_width),
                         &(renderer->_framebuffer_height));
  renderer->on_aspect_change();
}

} // namespace water_simulator::renderer