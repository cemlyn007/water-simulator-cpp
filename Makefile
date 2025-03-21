cpu:
	bazel run  //water_simulator/bin:water_simulator \
	 --subcommands -c opt

# It would be nice to find a way to not have to add the LD_LIBRARY_PATH for libOpenCL for CPU,
#  at the moment my system auto detects the CUDA one that leads to a runtime exception.
sycl_cpu:
	LD_LIBRARY_PATH="/opt/intel/oneapi/2025.0/lib" bazel run  //water_simulator/bin:water_simulator \
	 --subcommands -c opt --//water_simulator/engine:sycl=true --//water_simulator/engine:device=cpu

gpu:
	bazel run  //water_simulator/bin:water_simulator --subcommands -c opt \
	--//water_simulator/engine:sycl=true --//water_simulator/engine:device=nvidia

.PHONY: cpu sycl_cpu gpu
