# Water Simulator C++ (and SYCL*)

<p align="center">
  <img src="https://github.com/user-attachments/assets/4670c39a-85a4-4d30-8502-10061dae9fde">
</p>

This repository is basically a rewrite in C++ of my [Python Water Simulator](https://github.com/cemlyn007/water-simulator/tree/main) that used JAX and OpenGL.
This implementation tries to stay faithful to the Python one, mainly because I am curious to see how much I can beat it by!
This work is building upon the [SYCL toolchain](https://github.com/cemlyn007/rules_sycl) that I wrote.

So far, I have made the decision to use SYCL instead of CUDA for the increased flexibility of being able to easily implement parallelisation on the CPU or NVIDIA,
although in future work I may try a CUDA implementation to see if I find any performance benefits.
The CUDA implementation definitely has the appeal for interopability with OpenGL which I am yet to explore.

## Intended Goals

1. Increased proficiency in writing C++.
2. Use SYCL to reduce latency by offloading the simulator logic to the GPU.
3. Increased proficiency in writing Bazel.
4. Optimize performance by further familarising with using perf and nsys.

## Achieved Goals so far
1. The CPU non-SYCL (at least) has been confirmed to beat my original Python implementation for the **default** grid size, achieving locally approximately 600 microseconds per frame.

*My SYCL implementation is still in progress
