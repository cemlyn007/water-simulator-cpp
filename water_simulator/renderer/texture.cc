#include "water_simulator/renderer/texture.h"

namespace water_simulator::renderer {

Texture::Texture(GLuint texture) : _texture(texture) {}

Texture::~Texture() {}

void Texture::use() {
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, _texture);
}

void Texture::unuse() { glBindTexture(GL_TEXTURE_2D, 0); }

} // namespace water_simulator::renderer