#include "water_simulator/renderer/camera.h"
#include <stdexcept>

namespace water_simulator::renderer {
Camera::Camera(int width, int height) : _width(width), _height(height) {
  glGenFramebuffers(1, &_framebuffer);
  glBindFramebuffer(GL_FRAMEBUFFER, _framebuffer);

  // Create texture
  glGenTextures(1, &rendered_texture);
  glBindTexture(GL_TEXTURE_2D, rendered_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, _width, _height, 0, GL_RGB,
               GL_UNSIGNED_BYTE, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // Create depth buffer
  glGenRenderbuffers(1, &_depth_render_buffer);
  glBindRenderbuffer(GL_RENDERBUFFER, _depth_render_buffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, _width, _height);
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, _depth_render_buffer);

  // Attach texture
  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rendered_texture,
                       0);
  glDrawBuffer(GL_COLOR_ATTACHMENT0);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    throw std::runtime_error("Failed to create camera");
  }

  unbind();
}

Camera::~Camera() {
  glDeleteRenderbuffers(1, &_depth_render_buffer);
  glDeleteFramebuffers(1, &_framebuffer);
  glDeleteTextures(1, &rendered_texture);
}

void Camera::resize(int width, int height) {
  _width = width;
  _height = height;

  glBindFramebuffer(GL_FRAMEBUFFER, _framebuffer);

  // Resize texture
  glBindTexture(GL_TEXTURE_2D, rendered_texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, _width, _height, 0, GL_RGB,
               GL_UNSIGNED_BYTE, nullptr);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // Resize depth buffer
  glBindRenderbuffer(GL_RENDERBUFFER, _depth_render_buffer);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, _width, _height);

  // Reattach buffers
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, _depth_render_buffer);
  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, rendered_texture,
                       0);
}

void Camera::bind() {
  glBindFramebuffer(GL_FRAMEBUFFER, _framebuffer);
  glViewport(0, 0, _width, _height);
}

void Camera::unbind() { glBindFramebuffer(GL_FRAMEBUFFER, 0); }

Texture Camera::texture() { return Texture(rendered_texture); }

} // namespace water_simulator::renderer
