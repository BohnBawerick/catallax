set -uo pipefail
export CARGO_TERM_COLOR=always CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0 CARGO_PROFILE_TEST_OPT_LEVEL=2 SOLIDDRAFT_REQUIRE_GPU=1 LP_NUM_THREADS=1 WGPU_BACKEND=vulkan
echo "== hosted.yml partition-linux setup: mesa-vulkan-drivers"
SUDO=; [ "$(id -u)" -ne 0 ] && SUDO=sudo
$SUDO apt-get -o Acquire::Retries=2 -o Acquire::http::Timeout=30 update && $SUDO apt-get -o Acquire::Retries=2 -o Acquire::http::Timeout=30 install -y --no-install-recommends mesa-vulkan-drivers
echo "apt_exit=$?"
echo "== note: no tests/ribbon.rs target at 59f3fb2b; dropped --test ribbon as instructed"
status=0
for t in ux_options_bar ring_menu tool_palette; do
  echo "== cargo test -p soliddraft-app --test $t -- --test-threads=1"
  cargo test -p soliddraft-app --test $t -- --test-threads=1
  rc=$?; echo "TARGET_RESULT $t exit=$rc"; [ $rc -ne 0 ] && status=1
done
exit $status
