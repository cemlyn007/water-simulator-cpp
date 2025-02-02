#include "water_simulator/renderer/renderer.h"
#include <chrono>

using namespace std::chrono_literals;

int main(int argc, char *argv[]) {
  water_simulator::renderer::init();
  water_simulator::renderer::Renderer renderer(1080, 1080);
  while (!renderer.should_close()) {
    renderer.render();
  }
  water_simulator::renderer::terminate();
  return 0;
}