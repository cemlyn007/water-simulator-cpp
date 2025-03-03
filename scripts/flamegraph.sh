set -e

workspace=$(bazel info workspace)

tmp_directory=/tmp/water-simulator-cpp-flamegraph
perf_file_path=$tmp_directory/perf.data
svg_file_path=$tmp_directory/perf.svg

# -fno-omit-frame-pointer is advised from http://www.brendangregg.com/FlameGraphs/cpuflamegraphs.html#Instructions.
bazel build --compilation_mode=opt -c opt --copt="-g" --copt="-fno-omit-frame-pointer" //water_simulator/bin:water_simulator

mkdir -p $tmp_directory
echo "Require sudo to run perf:"
sudo perf record -a -b -g --output=$perf_file_path $workspace/bazel-bin/water_simulator/bin/water_simulator && sudo chown $USER $perf_file_path

perf script --input=$perf_file_path | third_party/FlameGraph/stackcollapse-perf.pl | $workspace/third_party/FlameGraph/flamegraph.pl > $svg_file_path

if [ -x "$(command -v google-chrome)" ]; then
  google-chrome $svg_file_path
fi
echo "Perf data recorded at $perf_file_path"
echo "Flamegraph generated at $svg_file_path"
