use core::cmp::Ordering;
use core::ops::{Add, Div, Mul, Sub};

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Orient {
    CounterClockwise,
    Clockwise,
    Collinear,
}

pub trait GlobalFunctions {
    type Point: Point<Number = Self::Number>;
    type Number: Number;

    fn dist2(&self, p: &Self::Point, q: &Self::Point) -> Self::Number;
    fn orient(&self, p: &Self::Point, q: &Self::Point, r: &Self::Point) -> Orient;
    fn circumdelta(
        &self,
        p: &Self::Point,
        b: &Self::Point,
        c: &Self::Point,
    ) -> (Self::Number, Self::Number);
    fn circumradius2(&self, p: &Self::Point, b: &Self::Point, c: &Self::Point) -> Self::Number;
    fn circumcenter(&self, p: &Self::Point, b: &Self::Point, c: &Self::Point) -> Self::Point;
    fn in_circle(&self, p: &Self::Point, b: &Self::Point, c: &Self::Point, q: &Self::Point)
        -> bool;
    fn nearly_equals(&self, p: &Self::Point, q: &Self::Point) -> bool;
    /// Can be overridden to use stable sorting. By default, this runs `slice::sort_unstable_by`.
    fn sort_slice_by<T>(&self, s: &mut [T], f: impl FnMut(&T, &T) -> Ordering);
}

#[cfg(feature = "robust")]
fn into_robust_coord<N: Number + Into<f64>>(p: &impl Point<Number = N>) -> robust::Coord<N> {
    robust::Coord { x: p.x(), y: p.y() }
}

#[cfg(feature = "robust")]
#[derive(Debug)]
pub struct DefaultGlobalFunctions<P: Point>(core::marker::PhantomData<P>);
#[cfg(feature = "robust")]
impl<P: Point> DefaultGlobalFunctions<P> {
    pub const fn const_new() -> Self {
        Self(core::marker::PhantomData)
    }
}
#[cfg(feature = "robust")]
impl<P: Point> Copy for DefaultGlobalFunctions<P> {}
#[cfg(feature = "robust")]
impl<P: Point> Clone for DefaultGlobalFunctions<P> {
    fn clone(&self) -> Self {
        *self
    }
}
#[cfg(feature = "robust")]
impl<P: Point> Default for DefaultGlobalFunctions<P> {
    fn default() -> Self {
        Self(core::marker::PhantomData)
    }
}
#[cfg(feature = "robust")]
impl<P: Point> GlobalFunctions for DefaultGlobalFunctions<P>
where
    P::Number: Into<f64>,
{
    type Point = P;
    type Number = P::Number;

    fn dist2(&self, p: &Self::Point, q: &Self::Point) -> Self::Number {
        let dx = p.x() - q.x();
        let dy = p.y() - q.y();
        dx * dx + dy * dy
    }

    fn orient(&self, p: &Self::Point, q: &Self::Point, r: &Self::Point) -> Orient {
        // Returns a **negative** value if ```self```, ```q``` and ```r``` occur in counterclockwise order (```r``` is to the left of the directed line ```self``` --> ```q```)
        // Returns a **positive** value if they occur in clockwise order(```r``` is to the right of the directed line ```self``` --> ```q```)
        // Returns zero is they are collinear
        match robust::orient2d(
            into_robust_coord(p),
            into_robust_coord(q),
            into_robust_coord(r),
        )
        .partial_cmp(&0.0)
        {
            Some(Ordering::Less) => Orient::CounterClockwise,
            Some(Ordering::Equal) => Orient::Collinear,
            Some(Ordering::Greater) => Orient::Clockwise,
            None => panic!("orient2d returned NaN"),
        }
    }

    fn circumdelta(
        &self,
        p: &Self::Point,
        b: &Self::Point,
        c: &Self::Point,
    ) -> (Self::Number, Self::Number) {
        let dx = b.x() - p.x();
        let dy = b.y() - p.y();
        let ex = c.x() - p.x();
        let ey = c.y() - p.y();

        let bl = dx * dx + dy * dy;
        let cl = ex * ex + ey * ey;
        let d = Self::Number::ONE_HALF / (dx * ey - dy * ex);

        let x = (ey * bl - dy * cl) * d;
        let y = (dx * cl - ex * bl) * d;
        (x, y)
    }

    fn circumradius2(&self, p: &Self::Point, b: &Self::Point, c: &Self::Point) -> Self::Number {
        let (x, y) = self.circumdelta(p, b, c);
        x * x + y * y
    }

    fn circumcenter(&self, p: &Self::Point, b: &Self::Point, c: &Self::Point) -> Self::Point {
        let (x, y) = self.circumdelta(p, b, c);
        Self::Point::new_point(p.x() + x, p.y() + y)
    }

    fn in_circle(
        &self,
        s: &Self::Point,
        b: &Self::Point,
        c: &Self::Point,
        p: &Self::Point,
    ) -> bool {
        let dx = s.x() - p.x();
        let dy = s.y() - p.y();
        let ex = b.x() - p.x();
        let ey = b.y() - p.y();
        let fx = c.x() - p.x();
        let fy = c.y() - p.y();

        let ap = dx * dx + dy * dy;
        let bp = ex * ex + ey * ey;
        let cp = fx * fx + fy * fy;

        dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx)
            < Self::Number::ZERO
    }

    fn nearly_equals(&self, s: &Self::Point, p: &Self::Point) -> bool {
        (s.x() - p.x()).abs() <= Self::Number::EPSILON
            && (s.y() - p.y()).abs() <= Self::Number::EPSILON
    }

    fn sort_slice_by<T>(&self, s: &mut [T], f: impl FnMut(&T, &T) -> Ordering) {
        s.sort_unstable_by(f)
    }
}

pub trait Point: Clone {
    type Number: Number;

    fn new_point(x: Self::Number, y: Self::Number) -> Self
    where
        Self: Sized;
    fn x(&self) -> Self::Number;
    fn y(&self) -> Self::Number;

    fn from_point(p: &impl Point<Number = Self::Number>) -> Self {
        Self::new_point(p.x(), p.y())
    }
    fn to_point<P: Point<Number = Self::Number>>(&self) -> P {
        P::from_point(self)
    }
    fn into_point<P: Point<Number = Self::Number>>(self) -> P {
        P::from_point(&self)
    }
}

pub trait Number:
    Copy
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

    const ZERO: Self;
    const ONE_HALF: Self;
    const ONE: Self;
    const TWO: Self;
    const THREE: Self;
    const FOUR: Self;

    const INFINITY: Self;
    const NEG_INFINITY: Self;

    fn abs(self) -> Self;
    fn floor(self) -> Self;
    fn sqrt(self) -> Self;

    fn from_usize_truncate(n: usize) -> Self;
    fn into_usize_truncate(self) -> usize;

    fn max(self, other: Self) -> Self;
    fn min(self, other: Self) -> Self;
}

#[cfg(feature = "robust")]
impl<N> Point for robust::Coord<N>
where
    N: Number + Into<f64>,
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
    fn from_usize_truncate(n: usize) -> Self {
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
    fn from_usize_truncate(n: usize) -> Self {
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
