cargo fmt --all -- --check && CARGO_BUILD_JOBS=4 cargo test --release -p catallax-measure && CARGO_BUILD_JOBS=4 cargo test --release -p catallax-app --test rational_lofts
