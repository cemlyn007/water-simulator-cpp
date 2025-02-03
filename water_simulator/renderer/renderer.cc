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
      _mouse_click(false), _escape_pressed(false), _light(),
      _container((100.0 * 0.02) / 2.0, 0.02 * 2.0) {

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);

  on_aspect_change();

  std::array<float, 3> light_position = {1.2, 4.0, 2.0};

  _light.set_color({1.0, 1.0, 1.0});
  _light.set_model(translate(scale(eye4d(), {0.2, 0.2, 0.2}), light_position));

  _container.set_color({0.7, 0.7, 0.7});
  _container.set_model(eye4d());

  _container.set_light_color({1.0, 1.0, 1.0});
  _container.set_light_position(light_position);

  _camera_position = {2.5, 3.535534, 2.5};
  _camera_radians[0] = 0.7853982;
  _camera_radians[1] = 0.7853982;

  update_camera();
}

Renderer::~Renderer() { glfwDestroyWindow(_window); };

void Renderer::render() {
  glfwMakeContextCurrent(_window);
  GL_CALL(glViewport(0, 0, _framebuffer_width, _framebuffer_height));
  GL_CALL(glClearColor(0.1, 0.1, 0.1, 1.0));
  GL_CALL(glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT));
  update_camera();
  _shader.use();
  _light.draw();
  _container.draw();
  GL_CALL(glfwSwapBuffers(_window));
  _last_mouse_position_in_pixels[0] = _mouse_position_in_pixels[0];
  _last_mouse_position_in_pixels[1] = _mouse_position_in_pixels[1];
  _scroll_offset = 0.0;
  GL_CALL(glfwPollEvents());
  _mouse_position_change_in_pixels[0] =
      _mouse_position_in_pixels[0] - _last_mouse_position_in_pixels[0];
  _mouse_position_change_in_pixels[1] =
      _mouse_position_in_pixels[1] - _last_mouse_position_in_pixels[1];
}

void Renderer::on_aspect_change() {
  float aspect = static_cast<float>(_framebuffer_width) /
                 static_cast<float>(_framebuffer_height);
  auto projection = perspective(radians(60), aspect, 0.01, 100.0);
  _light.set_projection(projection);
  _container.set_projection(projection);
};

void Renderer::update_camera() {
  float camera_radius =
      std::min(std::max(0.0, norm(_camera_position) + _scroll_offset), 25.0);
  _camera_radians[0] = std::fmod(
      _camera_radians[0] + radians(_mouse_position_change_in_pixels[0]),
      (2 * M_PI));
  _camera_radians[1] = std::fmod(
      _camera_radians[1] + radians(_mouse_position_change_in_pixels[1]),
      (2 * M_PI));
  _camera_position = update_orbit_camera_position(
      _camera_radians[0], _camera_radians[1], camera_radius);
  auto view = look_at(_camera_position, {0.0, 0.5, 0.0}, {0.0, 1.0, 0.0});
  _light.set_view(view);
  _container.set_view(view);
  _container.set_view_position(_camera_position);
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
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  glfwSetWindowUserPointer(window, this);
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GLFW_TRUE);
  glfwSetInputMode(window, GLFW_STICKY_MOUSE_BUTTONS, GLFW_TRUE);
  glfwSetKeyCallback(window, key_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetCursorPosCallback(window, cursor_position_callback);
  glfwSetScrollCallback(window, scroll_callback);
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

void Renderer::scroll_callback(GLFWwindow *window, double xoffset,
                               double yoffset) {
  Renderer *renderer =
      static_cast<Renderer *>(glfwGetWindowUserPointer(window));
  (renderer->_scroll_offset) = yoffset;
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