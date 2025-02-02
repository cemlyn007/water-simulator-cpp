#include "water_simulator/renderer/algebra.h"
#include <iostream>
#include <math.h>
#include <mdspan>

namespace water_simulator::renderer {

void print(const std::array<float, 16> &matrix) {
  for (const auto number : matrix) {
    std::cout << number << ", ";
  }
  std::cout << std::endl;
}

void print(const std::array<std::array<float, 4>, 4> &matrix) {
  for (const auto &row : matrix) {
    for (const auto number : row) {
      std::cout << number << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

std::array<float, 16> eye4d() {
  std::array<float, 16> result{};
  auto R_md = std::mdspan(result.data(), 4, 4);
  for (std::size_t i = 0; i < 4; ++i) {
    R_md[i, i] = 1.0;
  }
  return result;
}

std::array<float, 16> multiply_matrices(const std::array<float, 16> &a,
                                        const std::array<float, 16> &b) {
  std::array<float, 16> result;

  auto A_md = std::mdspan(a.data(), 4, 4);
  auto B_md = std::mdspan(b.data(), 4, 4);
  auto R_md = std::mdspan(result.data(), 4, 4);

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      float sum = 0.0f;
      for (std::size_t k = 0; k < 4; ++k) {
        sum += A_md[i, k] * B_md[k, j];
      }
      R_md[i, j] = sum;
    }
  }

  return result;
}

std::array<float, 16> translate(const std::array<float, 16> &matrix,
                                const std::array<float, 3> &vector) {
  // Create a 4x4 identity matrix for scaling.
  std::array<float, 16> translation = eye4d();

  auto translation_md = std::mdspan(translation.data(), 4, 4);

  // Replace the diagonal components with the scaling factors.
  translation_md[3, 0] = vector[0];
  translation_md[3, 1] = vector[1];
  translation_md[3, 2] = vector[2];

  // Multiply the input matrix by the scaling matrix.
  // This is equivalent to: result = matrix dot scaling_matrix.
  return multiply_matrices(matrix, translation);
}

std::array<float, 16> scale(const std::array<float, 16> &matrix,
                            const std::array<float, 3> &vector) {
  // Create a 4x4 identity matrix for scaling.
  std::array<float, 16> scaling = eye4d();

  auto scaling_md = std::mdspan(scaling.data(), 4, 4);

  // Replace the diagonal components with the scaling factors.
  scaling_md[0, 0] = vector[0];
  scaling_md[1, 1] = vector[1];
  scaling_md[2, 2] = vector[2];

  // Multiply the input matrix by the scaling matrix.
  // This is equivalent to: result = matrix dot scaling_matrix.
  return multiply_matrices(matrix, scaling);
}

std::array<float, 16> look_at(const std::array<float, 3> &eye,
                              const std::array<float, 3> &center,
                              const std::array<float, 3> &up) {
  // Compute forward vector f = normalize(center - eye)
  std::array<float, 3> f = {center[0] - eye[0], center[1] - eye[1],
                            center[2] - eye[2]};
  float f_norm = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
  f[0] /= f_norm;
  f[1] /= f_norm;
  f[2] /= f_norm;

  // Normalize the provided up vector to get u
  std::array<float, 3> u = up;
  float u_norm = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
  u[0] /= u_norm;
  u[1] /= u_norm;
  u[2] /= u_norm;

  // Compute side vector s = normalize(cross(f, u))
  std::array<float, 3> s = {f[1] * u[2] - f[2] * u[1],
                            f[2] * u[0] - f[0] * u[2],
                            f[0] * u[1] - f[1] * u[0]};
  float s_norm = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
  s[0] /= s_norm;
  s[1] /= s_norm;
  s[2] /= s_norm;

  u = {s[1] * f[2] - s[2] * f[1], s[2] * f[0] - s[0] * f[2],
       s[0] * f[1] - s[1] * f[0]};

  std::array<float, 16> flat_result = eye4d();

  auto result = std::mdspan(flat_result.data(), 4, 4);

  result[0, 0] = s[0];
  result[1, 0] = s[1];
  result[2, 0] = s[2];

  result[0, 1] = u[0];
  result[1, 1] = u[1];
  result[2, 1] = u[2];

  result[0, 2] = -f[0];
  result[1, 2] = -f[1];
  result[2, 2] = -f[2];

  result[3, 0] = -(s[0] * eye[0] + s[1] * eye[1] + s[2] * eye[2]);
  result[3, 1] = -(u[0] * eye[0] + u[1] * eye[1] + u[2] * eye[2]);
  result[3, 2] = f[0] * eye[0] + f[1] * eye[1] + f[2] * eye[2];

  return flat_result;
}

std::array<float, 16> perspective(float fov, float aspect, float near,
                                  float far) {
  float f = 1.0f / std::tan(fov / 2.0f);
  float nf = 1.0f / (near - far);

  std::array<float, 16> flat_result{};
  auto result = std::mdspan(flat_result.data(), 4, 4);

  result[0, 0] = f / aspect;
  result[1, 1] = f;
  result[2, 2] = (far + near) * nf;
  result[2, 3] = -1.0f;
  result[3, 2] = (2.0f * far * near) * nf;

  return flat_result;
}

} // namespace water_simulator::renderer
