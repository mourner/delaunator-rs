use core::cmp::Ordering;
use core::ops::{Add, Div, Mul, Sub};

/// Trait for types that are 2D points of a given number.
pub trait Point: Clone {
    type Number: Number;

    /// Creates a new point with the specified numbers.
    fn new_point(x: Self::Number, y: Self::Number) -> Self
    where
        Self: Sized;
    /// Returns the `x` coordinate.
    fn x(&self) -> Self::Number;
    /// Returns the `y` coordiante.
    fn y(&self) -> Self::Number;

    /// Creates a new point from another point type.
    #[inline]
    fn from_point(p: &impl Point<Number = Self::Number>) -> Self {
        Self::new_point(p.x(), p.y())
    }
    /// Converts `&self` into another point type.
    #[inline]
    fn to_point<P: Point<Number = Self::Number>>(&self) -> P {
        P::from_point(self)
    }
    /// Converts `self` into another point type.
    #[inline]
    fn into_point<P: Point<Number = Self::Number>>(self) -> P {
        P::from_point(&self)
    }

    /// Can be overridden to use stable sorting. By default, this runs `slice::sort_unstable_by`.
    #[inline]
    fn sort_slice_by<T>(s: &mut [T], f: impl FnMut(&T, &T) -> Ordering) {
        s.sort_unstable_by(f)
    }
}

/// A number that can be used in a triangulation point.
///
/// Implemented by default for `f64` and `f32`.
pub trait Number:
    Copy
    + Into<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + PartialOrd
{
    /// Near-duplicate points (where both `x` and `y` only differ within this value)
    /// will not be included in the triangulation for robustness.
    ///
    /// `std` implementations use `f{32,64}::EPSILON * 2.0`
    const EPSILON: Self;

    /// The constant value of `0.`.
    const ZERO: Self;
    /// The constant value of `0.5`.
    const ONE_HALF: Self;
    /// The constant value of `1.`.
    const ONE: Self;
    /// The constant value of `2.`.
    const TWO: Self;
    /// The constant value of `3.`.
    const THREE: Self;
    /// The constant value of `4.`.
    const FOUR: Self;

    /// Infinity
    ///
    /// std implementations use `f{32,64}::INFINITY`
    const INFINITY: Self;
    /// Negative Infinity
    ///
    /// std implementations use `f{32,64}::NEG_INFINITY`
    const NEG_INFINITY: Self;

    /// The absolute value of this number.
    fn abs(self) -> Self;
    /// The floor of this number.
    fn floor(self) -> Self;
    /// The squareroot of this number.
    fn sqrt(self) -> Self;

    /// Converts a `usize` into the nearest value.
    ///
    /// std implementations use `n as Self`
    fn from_usize(n: usize) -> Self;
    /// Converts `self` into a `usize` by truncating the decimal portion.
    ///
    /// std implementations use `self as usize`
    fn into_usize_truncate(self) -> usize;

    /// The maximum of this and another number.
    fn max(self, other: Self) -> Self;
    /// The minimum of this and another number.
    fn min(self, other: Self) -> Self;
}

impl<N> Point for robust::Coord<N>
where
    N: Number,
{
    type Number = N;

    #[inline]
    fn new_point(x: Self::Number, y: Self::Number) -> Self
    where
        Self: Sized,
    {
        Self { x, y }
    }

    #[inline]
    fn x(&self) -> Self::Number {
        self.x
    }

    #[inline]
    fn y(&self) -> Self::Number {
        self.y
    }
}
impl<N> Point for (N, N)
where
    N: Number,
{
    type Number = N;

    #[inline]
    fn new_point(x: Self::Number, y: Self::Number) -> Self {
        (x, y)
    }

    #[inline]
    fn x(&self) -> Self::Number {
        self.0
    }

    #[inline]
    fn y(&self) -> Self::Number {
        self.1
    }
}
impl<N> Point for [N; 2]
where
    N: Number,
{
    type Number = N;

    #[inline]
    fn new_point(x: Self::Number, y: Self::Number) -> Self {
        [x, y]
    }

    #[inline]
    fn x(&self) -> Self::Number {
        self[0]
    }

    #[inline]
    fn y(&self) -> Self::Number {
        self[1]
    }
}

impl Number for f64 {
    const EPSILON: Self = f64::EPSILON * 2.0;

    const ZERO: Self = 0.0;
    const ONE_HALF: Self = 0.5;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const THREE: Self = 3.0;
    const FOUR: Self = 4.0;

    const INFINITY: Self = f64::INFINITY;
    const NEG_INFINITY: Self = f64::NEG_INFINITY;

    #[cfg(feature = "std")]
    #[inline]
    fn abs(self) -> Self {
        f64::abs(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn abs(self) -> Self {
        const SIGN_BIT: u64 = 1 << 63;
        f64::from_bits(f64::to_bits(self) & !SIGN_BIT)
    }

    #[cfg(feature = "std")]
    #[inline]
    fn floor(self) -> Self {
        f64::floor(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn floor(self) -> Self {
        let mut res = (self as i64) as f64;
        if res > self {
            res -= 1.0;
        }
        res
    }

    #[cfg(feature = "std")]
    #[inline]
    fn sqrt(self) -> Self {
        f64::sqrt(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn sqrt(self) -> Self {
        if self < 2.0 {
            return self;
        };

        let sc = Number::sqrt(self / 4.0) * 2.0;
        let lc = sc + 1.0;

        if lc * lc > self {
            sc
        } else {
            lc
        }
    }

    #[inline]
    fn from_usize(n: usize) -> Self {
        n as Self
    }

    #[inline]
    fn into_usize_truncate(self) -> usize {
        self as usize
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f64::max(self, other)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f64::min(self, other)
    }
}
impl Number for f32 {
    const EPSILON: Self = f32::EPSILON * 2.0;

    const ZERO: Self = 0.0;
    const ONE_HALF: Self = 0.5;
    const ONE: Self = 1.0;
    const TWO: Self = 2.0;
    const THREE: Self = 3.0;
    const FOUR: Self = 4.0;

    const INFINITY: Self = f32::INFINITY;
    const NEG_INFINITY: Self = f32::NEG_INFINITY;

    #[cfg(feature = "std")]
    #[inline]
    fn abs(self) -> Self {
        f32::abs(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn abs(self) -> Self {
        const SIGN_BIT: u32 = 1 << 31;
        f32::from_bits(f32::to_bits(self) & !SIGN_BIT)
    }

    #[cfg(feature = "std")]
    #[inline]
    fn floor(self) -> Self {
        f32::floor(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn floor(self) -> Self {
        let mut res = (self as i32) as f32;
        if res > self {
            res -= 1.0;
        }
        res
    }

    #[cfg(feature = "std")]
    #[inline]
    fn sqrt(self) -> Self {
        f32::sqrt(self)
    }

    #[cfg(not(feature = "std"))]
    #[inline]
    fn sqrt(self) -> Self {
        if self < 2.0 {
            return self;
        };

        let sc = Number::sqrt(self / 4.0) * 2.0;
        let lc = sc + 1.0;

        if lc * lc > self {
            sc
        } else {
            lc
        }
    }

    #[inline]
    fn from_usize(n: usize) -> Self {
        n as Self
    }

    #[inline]
    fn into_usize_truncate(self) -> usize {
        self as usize
    }

    #[inline]
    fn max(self, other: Self) -> Self {
        f32::max(self, other)
    }

    #[inline]
    fn min(self, other: Self) -> Self {
        f32::min(self, other)
    }
}
