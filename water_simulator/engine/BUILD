load("@bazel_skylib//lib:selects.bzl", "selects")
load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//water_simulator:__subpackages__"])

cc_library(
    name = "engine",
    srcs = selects.with_or({
        ("//water_simulator:sycl_and_sycl_cpu", "//water_simulator:sycl_and_sycl_nvidia"): [
            "raycast.cc",
            "sycl_engine.cc",
        ],
        "//conditions:default": [
            "engine.cc",
            "raycast.cc",
        ],
    }),
    hdrs = ["engine.h"],
    cxxopts = select({
        "//water_simulator:sycl_and_sycl_nvidia": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64,nvptx64-nvidia-cuda",
        ],
        "//water_simulator:sycl_and_sycl_cpu": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64",
        ],
        "//conditions:default": [
        ],
    }),
    linkopts = select({
        "//water_simulator:sycl_and_sycl_nvidia": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64,nvptx64-nvidia-cuda",
            "--offload-arch=sm_80",
        ],
        "//water_simulator:sycl_and_sycl_cpu": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64",
        ],
        "//conditions:default": [
        ],
    }),
    deps = [":state"] + selects.with_or({
        ("//water_simulator:sycl_and_sycl_cpu", "//water_simulator:sycl_and_sycl_nvidia"): [
            "@local_config_sycl//sycl",
        ],
        "//conditions:default": [
        ],
    }),
)

cc_library(
    name = "state",
    srcs = ["state.cc"],
    hdrs = ["state.h"],
)
