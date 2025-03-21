load("@bazel_skylib//lib:selects.bzl", "selects")
load("@bazel_skylib//rules:common_settings.bzl", "bool_flag", "string_flag")
load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//water_simulator:__subpackages__"])

bool_flag(
    name = "sycl",
    build_setting_default = False,
    visibility = ["//visibility:public"],
)

config_setting(
    name = "sycl_enabled",
    flag_values = {
        ":sycl": "true",
    },
)

string_flag(
    name = "device",
    build_setting_default = "cpu",
    visibility = ["//visibility:public"],
)

config_setting(
    name = "cpu",
    flag_values = {
        ":device": "cpu",
    },
)

config_setting(
    name = "nvidia",
    flag_values = {
        ":device": "nvidia",
    },
)

selects.config_setting_group(
    name = "sycl_and_sycl_nvidia",
    match_all = [
        ":sycl_enabled",
        ":nvidia",
    ],
)

selects.config_setting_group(
    name = "sycl_and_sycl_cpu",
    match_all = [
        ":sycl_enabled",
        ":cpu",
    ],
)

cc_library(
    name = "engine",
    srcs = selects.with_or({
        (":sycl_and_sycl_cpu", "sycl_and_sycl_nvidia"): [
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
        ":sycl_and_sycl_cpu": [
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
        ":sycl_and_sycl_cpu": [
            "-fsycl",
            "-fsycl-unnamed-lambda",
            "-fsycl-targets=spir64",
        ],
        "//conditions:default": [
        ],
    }),
    deps = [":state"] + selects.with_or({
        (":sycl_and_sycl_cpu", "sycl_and_sycl_nvidia"): [
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
