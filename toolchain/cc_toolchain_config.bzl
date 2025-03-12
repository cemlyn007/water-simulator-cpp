load("@bazel_tools//tools/build_defs/cc:action_names.bzl", "ACTION_NAMES")
load(
    "@bazel_tools//tools/cpp:cc_toolchain_config_lib.bzl",
    "action_config",
    "artifact_name_pattern",
    "env_entry",
    "env_set",
    "feature",
    "feature_set",
    "flag_group",
    "flag_set",
    "tool",
    "tool_path",
    "variable_with_value",
    "with_feature_set",
)


def _impl(ctx):
    tool_paths = [
        tool_path(
            name = "gcc",
            path = "/opt/intel/oneapi/2025.0/bin/icpx",
        ),
        tool_path(
            name = "ld",
            path = "/usr/bin/ld",
        ),
        tool_path(
            name = "ar",
            path = "/usr/bin/ar",
        ),
        tool_path(
            name = "cpp",
            path = "/bin/false",
        ),
        tool_path(
            name = "gcov",
            path = "/bin/false",
        ),
        tool_path(
            name = "nm",
            path = "/bin/false",
        ),
        tool_path(
            name = "objdump",
            path = "/bin/false",
        ),
        tool_path(
            name = "strip",
            path = "/bin/false",
        ),
    ]

    action_configs = [
        action_config(
            ACTION_NAMES.c_compile,
            tools = [tool(path = "/usr/bin/gcc")],
        ),
        action_config(
            ACTION_NAMES.cpp_compile,
            tools = [tool(path = "/opt/intel/oneapi/2025.0/bin/icpx")],
        ),
    ]

    return cc_common.create_cc_toolchain_config_info(
        ctx = ctx,
        action_configs=action_configs,
        cxx_builtin_include_directories = [
            "/opt/intel/oneapi/2025.0/include",
            "/opt/intel/oneapi/compiler/2025.0/include",
            "/usr/include",
            "/opt/intel/oneapi/compiler/2025.0/lib/clang/19/include",
            "/opt/intel/oneapi/compiler/2025.0/opt/compiler/include",
            "/usr/lib/gcc/x86_64-linux-gnu/14/include"
        ],
        toolchain_identifier = "local",
        host_system_name = "local",
        target_system_name = "local",
        target_cpu = "k8",
        target_libc = "unknown",
        compiler = "clang",
        abi_version = "unknown",
        abi_libc_version = "unknown",
        tool_paths = tool_paths,
    )

cc_toolchain_config = rule(
    implementation = _impl,
    attrs = {},
    provides = [CcToolchainConfigInfo],
)
