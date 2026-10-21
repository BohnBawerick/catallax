//! Arithmetic forwarding to each type's exact implementation.

macro_rules! binops {
    ($t:ty; $($tr:ident :: $m:ident),+ $(,)?) => { $(
        impl ::std::ops::$tr for $t {
            type Output = $t;
            fn $m(self, rhs: $t) -> $t { <$t>::$m(&self, &rhs) }
        }
        impl ::std::ops::$tr<&$t> for &$t {
            type Output = $t;
            fn $m(self, rhs: &$t) -> $t { <$t>::$m(self, rhs) }
        }
    )+ };
}
macro_rules! negop {
    ($t:ty) => {
        impl ::std::ops::Neg for $t {
            type Output = $t;
            fn neg(self) -> $t {
                <$t>::neg(&self)
            }
        }
        impl ::std::ops::Neg for &$t {
            type Output = $t;
            fn neg(self) -> $t {
                <$t>::neg(self)
            }
        }
    };
}
