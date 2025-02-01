#ifndef WATER_SIMULATOR_RENDERER_SHADER_H_
#define WATER_SIMULATOR_RENDERER_SHADER_H_
#include <GL/glew.h>
#include <stdlib.h>
#include <string>

namespace water_simulator::renderer {

class Shader {
private:
  GLuint _program;

public:
  Shader(std::string vertex_source_code, std::string fragment_source_code);
  ~Shader();
  void use();
  void set_uniform(std::string name, float x, float y);

private:
  GLuint load_program(std::string vertex_source_code,
                      std::string fragment_source_code);
  GLuint load_vertex_shader(std::string vertex_source_code);
  GLuint load_fragment_shader(std::string fragment_source);
};

std::string read_file(const std::string &file_path);

} // namespace water_simulator::renderer
#endif