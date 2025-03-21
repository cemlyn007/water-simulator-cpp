#include <vector>
namespace water_simulator::renderer::entities {
void update_water_normals(std::vector<float> &vertex_normals, std::vector<float> &face_normals,
                          const std::vector<float> &heights, size_t resolution, const std::vector<float> &xz,
                          const std::vector<unsigned int> &indices, const std::vector<size_t> &count);

}