load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//water_simulator:__subpackages__"])

cc_library(
    name = "engine",
    srcs = select({
        "//conditions:default": [
            "engine.cc",
            "raycast.cc",
        ],
    }),
    hdrs = ["engine.h"],
    deps = [":state"],
)

cc_library(
    name = "state",
    srcs = ["state.cc"],
    hdrs = ["state.h"],
)
