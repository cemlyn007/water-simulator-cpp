#pragma once
#include <GL/glew.h>

namespace water_simulator::renderer {

class Texture {
private:
  GLuint _texture;

public:
  Texture(GLuint _texture);
  ~Texture();

  void use();
  void unuse();
};

} // namespace water_simulator::renderer