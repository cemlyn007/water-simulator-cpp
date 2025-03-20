load("@bazel_skylib//lib:selects.bzl", "selects")
load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//water_simulator:__subpackages__"])

config_setting(
    name = "sycl",
    values = {
        "define": "sycl_engine=true",
    },
)

config_setting(
    name = "sycl_nvidia",
    values = {
        "define": "sycl_nvidia=true",
    },
)

selects.config_setting_group(
    name = "sycl_and_sycl_nvidia",
    match_all = [
        ":sycl",
        ":sycl_nvidia",
    ],
)

cc_library(
    name = "engine",
    srcs = select({
        ":sycl": [
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
        ":sycl_and_sycl_nvidia": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64,nvptx64-nvidia-cuda",
        ],
        ":sycl": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64",
        ],
        "//conditions:default": [
        ],
    }),
    linkopts = select({
        ":sycl_and_sycl_nvidia": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64,nvptx64-nvidia-cuda",
            "--offload-arch=sm_80",
        ],
        ":sycl": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64",
        ],
        "//conditions:default": [
        ],
    }),
    deps = [":state"] + select({
        ":sycl": [
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
