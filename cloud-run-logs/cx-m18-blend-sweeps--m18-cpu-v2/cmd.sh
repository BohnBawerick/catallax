cargo fmt --all --check && cargo clippy -p catallax-field -p catallax-compare -p catallax-cli --all-targets -- -D warnings && cargo test --workspace --release --exclude catallax-preview
