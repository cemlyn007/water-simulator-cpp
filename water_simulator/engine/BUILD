load("@bazel_skylib//lib:selects.bzl", "selects")
load("@rules_cc//cc:defs.bzl", "cc_library")
load("//water_simulator:variables.bzl", "SYCL_CXXOPTS", "SYCL_DEPS", "SYCL_LINKOPTS")

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
    cxxopts = SYCL_CXXOPTS,
    linkopts = SYCL_LINKOPTS,
    deps = [":state"] + SYCL_DEPS,
)

cc_library(
    name = "state",
    srcs = ["state.cc"],
    hdrs = ["state.h"],
)
