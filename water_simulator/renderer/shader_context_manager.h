#pragma once
#include "water_simulator/renderer/shader.h"

namespace water_simulator::renderer {

class ShaderContextManager {
public:
  ShaderContextManager(Shader &shader);
  ~ShaderContextManager();

private:
  Shader &_shader;
};

} // namespace water_simulator::renderer