"""Variables for the water_simulator package."""

load("@bazel_skylib//lib:selects.bzl", "selects")

SYCL_CXXOPTS = select({
    "//water_simulator:sycl_and_sycl_nvidia": [
        "-fsycl",
        "-fsycl-unnamed-lambda",
        "-fsycl-targets=nvptx64-nvidia-cuda",
        "-Xsycl-target-backend",
        "--cuda-gpu-arch=sm_89",
    ],
    "//water_simulator:sycl_and_sycl_cpu": [
        "-fsycl",
        "-fsycl-unnamed-lambda",
        "-fsycl-targets=spir64",
    ],
    "//conditions:default": [
    ],
})

SYCL_LINKOPTS = SYCL_CXXOPTS + select({
    "//water_simulator:sycl_and_sycl_nvidia": [
    ],
    "//conditions:default": [
    ],
})

SYCL_DEPS = selects.with_or({
    ("//water_simulator:sycl_and_sycl_cpu", "//water_simulator:sycl_and_sycl_nvidia"): [
        "@local_config_sycl//sycl",
    ],
    "//conditions:default": [
    ],
})
