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
    name = "renderer",
    srcs = ["renderer.cc"],
    hdrs = ["renderer.h"],
    data = [
        "//water_simulator/renderer/shaders:shader.frag",
        "//water_simulator/renderer/shaders:shader.vert",
    ],
    linkopts = opengl_linkopts,
    deps = [
        ":gl_error_macro",
        ":shader",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "shader",
    srcs = ["shader.cc"],
    hdrs = ["shader.h"],
    linkopts = opengl_linkopts,
    deps = [
        ":gl_error_macro",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "gl_error_macro",
    srcs = ["gl_error_macro.cc"],
    hdrs = ["gl_error_macro.h"],
    linkopts = opengl_linkopts,
    deps = [
        "@glew//:glew_static",
        "@glfw",
    ],
)
