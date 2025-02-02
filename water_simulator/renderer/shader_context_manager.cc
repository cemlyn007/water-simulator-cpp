#include "water_simulator/renderer/shader_context_manager.h"

namespace water_simulator::renderer {

ShaderContextManager::ShaderContextManager(Shader &shader) : _shader(shader) {
  _shader.use();
}

ShaderContextManager::~ShaderContextManager() { _shader.unuse(); }

} // namespace water_simulator::renderer
