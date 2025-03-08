#include "water_simulator/renderer/shader.h"
#include "water_simulator/renderer/gl_error_macro.h"
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace water_simulator::renderer {

Shader::Shader(std::string vertex_source_code, std::string fragment_source_code)
    : _program(load_program(vertex_source_code, fragment_source_code)) {}

GLuint Shader::load_program(std::string vertex_source_code, std::string fragment_source_code) {
  GLuint vertex_shader = load_vertex_shader(vertex_source_code);
  GLuint fragment_shader = load_fragment_shader(fragment_source_code);
  GLuint shader_program = glCreateProgram();
  GL_CALL(glAttachShader(shader_program, vertex_shader));
  GL_CALL(glAttachShader(shader_program, fragment_shader));
  GL_CALL(glLinkProgram(shader_program));
  glDeleteShader(vertex_shader);
  glDeleteShader(fragment_shader);
  int success;
  char info_log[512];
  GL_CALL(glGetProgramiv(shader_program, GL_LINK_STATUS, &success));
  if (!success) {
    glGetProgramInfoLog(shader_program, 512, nullptr, info_log);
    throw std::runtime_error(std::string("ERROR::SHADER::PROGRAM::FAILED\n") + info_log);
  }
  // else...
  return shader_program;
};

Shader::~Shader() {
  if (_program != 0)
    glDeleteProgram(_program);
}

void Shader::use() { GL_CALL(glUseProgram(_program)); }

void Shader::unuse() { GL_CALL(glUseProgram(0)); }

void Shader::set_uniform(const std::string &name, int value) {
  GLuint location = glGetUniformLocation(_program, name.c_str());
  GL_CALL(glUniform1i(location, value););
}

void Shader::set_uniform_vector(const std::string &name, const std::array<float, 2> &vector) {
  GLuint location = glGetUniformLocation(_program, name.c_str());
  GL_CALL(glUniform2fv(location, 1, vector.data()););
}

void Shader::set_uniform_vector(const std::string &name, const std::array<float, 3> &vector) {
  GLuint location = glGetUniformLocation(_program, name.c_str());
  GL_CALL(glUniform3fv(location, 1, vector.data()););
}

void Shader::set_uniform_matrix(const std::string &name, const std::array<float, 16> &matrix) {
  GLuint location = glGetUniformLocation(_program, name.c_str());
  GL_CALL(glUniformMatrix4fv(location, 1, GL_TRUE, matrix.data()););
}

GLuint Shader::load_vertex_shader(std::string vertex_source_code) {
  const char *vertex_shader_c_str = vertex_source_code.c_str();
  GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
  GL_CALL(glShaderSource(vertex_shader, 1, &vertex_shader_c_str, nullptr));
  GL_CALL(glCompileShader(vertex_shader));
  int success;
  char info_log[512];
  GL_CALL(glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success));
  if (!success) {
    glGetShaderInfoLog(vertex_shader, 512, nullptr, info_log);
    throw std::runtime_error(std::string("ERROR::SHADER::VERTEX::COMPILATION_FAILED\n") + info_log);
  }
  // else...
  return vertex_shader;
};

GLuint Shader::load_fragment_shader(std::string fragment_source) {
  const char *fragment_shader_source_c_str = fragment_source.c_str();
  GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
  GL_CALL(glShaderSource(fragment_shader, 1, &fragment_shader_source_c_str, nullptr));
  GL_CALL(glCompileShader(fragment_shader));
  int success;
  char info_log[512];
  GL_CALL(glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success));
  if (!success) {
    glGetShaderInfoLog(fragment_shader, 512, nullptr, info_log);
    throw std::runtime_error(std::string("ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n") + info_log);
  }
  // else...
  return fragment_shader;
};

std::string read_file(const std::string &file_path) {
  std::ifstream file(file_path);
  if (!file) {
    throw std::runtime_error("Could not open file: " + file_path);
  }
  // else...
  std::stringstream buffer;
  buffer << file.rdbuf();
  return buffer.str();
}

} // namespace water_simulator::renderer