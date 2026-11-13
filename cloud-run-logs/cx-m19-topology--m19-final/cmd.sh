cargo test --release -p catallax-cli --test m19_print_verdict -j 2 -- --nocapture --test-threads=1 && cargo fmt --all --check && git diff --check
