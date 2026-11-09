use catallax_numbers::Rational;

/// Classify every legal decimal rounding, independently of the STEP writer.
pub fn line_floor_bounds(lines: &[[[Rational; 3]; 2]], digits: u32) -> (bool, bool) {
    let floor = Rational::parse("1e-14").unwrap();
    let mut must_refuse = false;
    let mut may_refuse = false;
    for line in lines {
        let bounds = line.each_ref().map(|point| {
            point.each_ref().map(|coordinate| {
                let (lo, hi) = coordinate.to_decimal_bounds(digits.max(1));
                [Rational::parse(&lo).unwrap(), Rational::parse(&hi).unwrap()]
            })
        });
        let mut minimum = Rational::zero();
        let mut maximum = Rational::zero();
        for axis in 0..3 {
            let [a, b] = bounds.each_ref().map(|point| &point[axis]);
            let low = (&a[0] - &b[1]).max(&b[0] - &a[1]).max(Rational::zero());
            let high = (&a[0] - &b[1]).abs().max((&a[1] - &b[0]).abs());
            minimum = minimum + &low * &low;
            maximum = maximum + &high * &high;
        }
        must_refuse |= maximum <= floor;
        may_refuse |= minimum <= floor;
    }
    (must_refuse, may_refuse)
}
