#include "water_simulator/renderer/algebra.h"
#include <array>
#include <math.h>
#include <mdspan>

namespace water_simulator::renderer {

float radians(float degrees) { return (degrees * M_PI) / 180.0; };

float norm(std::array<float, 3> vector) {
  return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

std::array<float, 3> normalize(const std::array<float, 3> &vector) {
  float denominator = norm(vector);
  return {
      vector[0] / denominator,
      vector[1] / denominator,
      vector[2] / denominator,
  };
};

std::array<float, 3> update_orbit_camera_position(float azimuth_radians, float elevation_radians, float radius) {
  return {radius * std::cos(azimuth_radians) * std::cos(elevation_radians), radius * std::sin(elevation_radians),
          radius * std::sin(azimuth_radians) * std::cos(elevation_radians)};
};

std::array<float, 16> eye4d() {
  std::array<float, 16> result;
  std::fill(result.begin(), result.end(), 0.0);
  auto R_md = std::mdspan(result.data(), 4, 4);
  for (std::size_t i = 0; i < 4; ++i) {
    R_md[i, i] = 1.0;
  }
  return result;
}

std::array<float, 16> multiply_matrices(const std::array<float, 16> &a, const std::array<float, 16> &b) {
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

std::array<float, 4> multiply_matrix(const std::array<float, 16> &a, const std::array<float, 4> &b) {
  std::array<float, 4> result;
  auto A_md = std::mdspan(a.data(), 4, 4);
  for (std::size_t i = 0; i < 4; ++i) {
    float sum = A_md[i, 0] * b[0];
    for (std::size_t j = 1; j < 4; ++j) {
      sum += A_md[i, j] * b[j];
    }
    result[i] = sum;
  }
  return result;
};

std::array<float, 16> translate(const std::array<float, 16> &matrix, const std::array<float, 3> &vector) {
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

std::array<float, 16> scale(const std::array<float, 16> &matrix, const std::array<float, 3> &vector) {
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

std::array<float, 16> transpose(const std::array<float, 16> &matrix) {
  std::array<float, 16> transposed{};

  // Iterate through all elements
  for (size_t row = 0; row < 4; ++row) {
    for (size_t col = 0; col < 4; ++col) {
      // Original index: row * 4 + col
      // Transposed index: col * 4 + row
      transposed[col * 4 + row] = matrix[row * 4 + col];
    }
  }

  return transposed;
}

std::array<float, 16> look_at(const std::array<float, 3> &eye, const std::array<float, 3> &center,
                              const std::array<float, 3> &up) {
  // Compute forward vector f = normalize(center - eye)
  std::array<float, 3> f = {center[0] - eye[0], center[1] - eye[1], center[2] - eye[2]};
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
  std::array<float, 3> s = {f[1] * u[2] - f[2] * u[1], f[2] * u[0] - f[0] * u[2], f[0] * u[1] - f[1] * u[0]};
  float s_norm = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
  s[0] /= s_norm;
  s[1] /= s_norm;
  s[2] /= s_norm;

  u = {s[1] * f[2] - s[2] * f[1], s[2] * f[0] - s[0] * f[2], s[0] * f[1] - s[1] * f[0]};

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

  return transpose(flat_result);
}

std::array<float, 16> perspective(float fov, float aspect, float near, float far) {
  float f = 1.0f / std::tan(fov / 2.0f);
  float nf = 1.0f / (near - far);

  std::array<float, 16> flat_result{};
  auto result = std::mdspan(flat_result.data(), 4, 4);

  result[0, 0] = f / aspect;
  result[1, 1] = f;
  result[2, 2] = (far + near) * nf;
  result[2, 3] = -1.0f;
  result[3, 2] = (2.0f * far * near) * nf;

  return transpose(flat_result);
}

std::array<float, 16> inverse(std::array<float, 16> matrix) {
  constexpr int n = 4;
  std::array<float, 16> inv = eye4d();
  constexpr float threshold = 1e-8;

  for (int i = 0; i < n; i++) {
    // Find the pivot row for column i.
    int pivot = i;
    float max_val = std::abs(matrix[i * n + i]);
    for (int j = i + 1; j < n; j++) {
      float val = std::abs(matrix[j * n + i]);
      if (val > max_val) {
        max_val = val;
        pivot = j;
      }
    }
    if (std::abs(matrix[pivot * n + i]) < threshold) {
      // Matrix is singular; return identity or handle error.
      return eye4d();
    }
    // Swap the current row with the pivot row if needed.
    if (pivot != i) {
      for (int j = 0; j < n; j++) {
        std::swap(matrix[i * n + j], matrix[pivot * n + j]);
        std::swap(inv[i * n + j], inv[pivot * n + j]);
      }
    }
    // Divide the pivot row by the pivot value.
    float diag = matrix[i * n + i];
    for (int j = 0; j < n; j++) {
      matrix[i * n + j] /= diag;
      inv[i * n + j] /= diag;
    }
    // Eliminate all other rows.
    for (int k = 0; k < n; k++) {
      if (k != i) {
        float factor = matrix[k * n + i];
        for (int j = 0; j < n; j++) {
          matrix[k * n + j] -= factor * matrix[i * n + j];
          inv[k * n + j] -= factor * inv[i * n + j];
        }
      }
    }
  }
  return inv;
}

} // namespace water_simulator::renderer
