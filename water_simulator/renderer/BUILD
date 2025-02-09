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
    linkopts = opengl_linkopts,
    deps = [
        ":algebra",
        ":gl_error_macro",
        ":shader",
        "//water_simulator/renderer:camera",
        "//water_simulator/renderer/entities:ball",
        "//water_simulator/renderer/entities:container",
        "//water_simulator/renderer/entities:light",
        "//water_simulator/renderer/entities:water",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "camera",
    srcs = ["camera.cc"],
    hdrs = ["camera.h"],
    linkopts = opengl_linkopts,
    deps = [
        ":gl_error_macro",
        ":texture",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "algebra",
    srcs = ["algebra.cc"],
    hdrs = ["algebra.h"],
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
    name = "shader_context_manager",
    srcs = ["shader_context_manager.cc"],
    hdrs = ["shader_context_manager.h"],
    linkopts = opengl_linkopts,
    deps = [
        ":gl_error_macro",
        ":shader",
        "@glew//:glew_static",
        "@glfw",
    ],
)

cc_library(
    name = "texture",
    srcs = ["texture.cc"],
    hdrs = ["texture.h"],
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
