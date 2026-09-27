#pragma once

// Fused vector operations for PDLP on GPU.
// Double precision throughout.

#include "device.hpp"

namespace samaya {
namespace gpu {

// Projection onto bounds: x = clamp(x, lower, upper)
// n elements
void project_bounds(int n, double* x, const double* lower, const double* upper,
                    cudaStream_t stream = 0);

// axpby: y = a*x + b*y
// n elements
void axpby(int n, double a, const double* x, double b, double* y, cudaStream_t stream = 0);

// axpy: y += a*x
// n elements
void axpy(int n, double a, const double* x, double* y, cudaStream_t stream = 0);

// Scale: x *= a
// n elements
void scale(int n, double a, double* x, cudaStream_t stream = 0);

// Copy: y = x
// n elements
void copy(int n, const double* x, double* y, cudaStream_t stream = 0);

// Dot product: returns x^T * y
// n elements
double dot(int n, const double* x, const double* y, cudaStream_t stream = 0);

// Norm: returns ||x||_2
// n elements
double norm2(int n, const double* x, cudaStream_t stream = 0);

// Norm squared: returns ||x||_2^2
// n elements
double norm2_squared(int n, const double* x, cudaStream_t stream = 0);

// Inf norm: returns max |x_i|
// n elements
double norm_inf(int n, const double* x, cudaStream_t stream = 0);

// Sum: returns sum(x_i)
// n elements
double sum(int n, const double* x, cudaStream_t stream = 0);

// Primal residual: ||A*x - r|| / (1 + ||A||*||x|| + ||r||)
// We compute the numerator on device: ||A*x - r||
// The denominator is computed on host
double primal_residual_norm(int m, const double* Ax, const double* r, cudaStream_t stream = 0);

// Dual residual: ||A^T*y + d - c|| / (1 + ||A||*||y|| + ||c||)
// Compute numerator: ||A^T*y + d - c||
double dual_residual_norm(int n, const double* ATy, const double* d, const double* c,
                          cudaStream_t stream = 0);

// Gap: |c^T*x - b^T*y| / (1 + |c^T*x| + |b^T*y|)
// Compute numerator components
double gap_numerator(int n, int m, const double* c, const double* x, const double* b,
                     const double* y, cudaStream_t stream = 0);

}  // namespace gpu
}  // namespace samaya