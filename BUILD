load("@hedron_compile_commands//:refresh_compile_commands.bzl", "refresh_compile_commands")

refresh_compile_commands(
    name = "refresh_compile_commands",
    targets = {
        "//water_simulator/bin:water_simulator": "",
    },
)

refresh_compile_commands(
    name = "refresh_compile_commands_sycl",
    targets = {
        "//water_simulator/bin:water_simulator": "-c opt --//water_simulator:sycl=true --//water_simulator:device=nvidia",
    },
)
