load("@bazel_skylib//lib:selects.bzl", "selects")
load("@bazel_skylib//rules:common_settings.bzl", "bool_flag", "string_flag")

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
