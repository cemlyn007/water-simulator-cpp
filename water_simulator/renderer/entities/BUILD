load("@rules_cc//cc:defs.bzl", "cc_library")

package(default_visibility = ["//water_simulator:__subpackages__"])

opengl_linkopts = select({
    "@platforms//os:windows": ["-DEFAULTLIB:opengl32.lib"],
    "@platforms//os:linux": [
        "-lGL",
        "-lX11",
    ],
    "@platforms//os:osx": [
        "-framework",
        "OpenGL",
    ],
})

cc_library(
    name = "ball",
    srcs = ["ball.cc"],
    hdrs = ["ball.h"],
    data = ["//water_simulator/renderer/shaders:simple"],
    linkopts = opengl_linkopts,
    deps = [
        "//water_simulator/renderer:gl_error_macro",
        "//water_simulator/renderer:shader",
        "//water_simulator/renderer:shader_context_manager",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "container",
    srcs = ["container.cc"],
    hdrs = ["container.h"],
    data = ["//water_simulator/renderer/shaders:simple"],
    linkopts = opengl_linkopts,
    deps = [
        "//water_simulator/renderer:gl_error_macro",
        "//water_simulator/renderer:shader",
        "//water_simulator/renderer:shader_context_manager",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "light",
    srcs = ["light.cc"],
    hdrs = ["light.h"],
    data = ["//water_simulator/renderer/shaders:light_cube"],
    linkopts = opengl_linkopts,
    deps = [
        "//water_simulator/renderer:gl_error_macro",
        "//water_simulator/renderer:shader",
        "//water_simulator/renderer:shader_context_manager",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "water",
    srcs = ["water.cc"],
    hdrs = ["water.h"],
    data = ["//water_simulator/renderer/shaders:basic_lighting"],
    linkopts = opengl_linkopts,
    deps = [
        "//water_simulator/renderer:algebra",
        "//water_simulator/renderer:gl_error_macro",
        "//water_simulator/renderer:shader",
        "//water_simulator/renderer:shader_context_manager",
        "//water_simulator/renderer:texture",
        "@glew//:glew_static",
        "@glfw",
    ],
)
