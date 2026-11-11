// Freeze the M25 constructor as portable exact-rational Rust input.
// Build with -std=c++17 -O0 -ffp-contract=off, then rustfmt the output.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

int main() {
    std::puts("//! Frozen binary64 evaluation of the M25 wavy control formula, read as exact rationals.");
    std::puts("use super::{Q, q};");
    std::puts("const HEIGHTS: [[[u64; 6]; 6]; 2] = [");
    for (int side = 0; side < 2; ++side) {
        const double amp = side ? -0.8 : 1.0, shift = side ? 1.1 : 0.3;
        std::puts("[");
        for (int i = 0; i < 6; ++i) {
            std::printf("[");
            for (int j = 0; j < 6; ++j) {
                const double x = i, y = j;
                const double z = amp * std::sin(0.9*x + shift)
                    * std::cos(0.7*y - shift) + 0.05*x;
                uint64_t bits;
                std::memcpy(&bits, &z, sizeof bits);
                std::printf("0x%016llx,", static_cast<unsigned long long>(bits));
            }
            std::puts("],");
        }
        std::puts("],");
    }
    std::puts(R"(];
pub(super) fn pair() -> [Vec<Vec<[Q; 3]>>; 2] {
    HEIGHTS.map(|net| net.iter().enumerate().map(|(i, row)| {
        row.iter().enumerate().map(|(j, bits)| [
            q(i as i64), q(j as i64), Q::from_f64_exact(f64::from_bits(*bits)).unwrap()
        ]).collect()
    }).collect())
})");
}
