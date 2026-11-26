set -e
nproc; free -g; df -h .
t0=$(date +%s)
cargo build -p catallax-cli --release
t1=$(date +%s); echo "BUILD_SECONDS $((t1-t0))"
cargo test -p catallax-core --release
t2=$(date +%s); echo "TEST_SECONDS $((t2-t1))"
