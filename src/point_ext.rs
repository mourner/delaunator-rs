use core::cmp::Ordering;

use crate::{Number, Point};

#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash)]
pub enum Orient {
    CounterClockwise,
    Clockwise,
    Collinear,
}

/// Internal trait for additional point functions.
pub trait PointExt: Point {
    fn dist2(&self, q: &Self) -> Self::Number;
    fn orient(&self, q: &Self, r: &Self) -> Orient;
    fn circumdelta(&self, b: &Self, c: &Self) -> (Self::Number, Self::Number);
    fn circumradius2(&self, b: &Self, c: &Self) -> Self::Number;
    fn circumcenter(&self, b: &Self, c: &Self) -> Self;
    fn in_circle(&self, b: &Self, c: &Self, q: &Self) -> bool;
    fn nearly_equals(&self, q: &Self) -> bool;
}

fn into_robust_coord<N: Number>(p: &impl Point<Number = N>) -> robust::Coord<N> {
    robust::Coord { x: p.x(), y: p.y() }
}

impl<P: Point> PointExt for P {
    fn dist2(&self, q: &Self) -> Self::Number {
        let dx = self.x() - q.x();
        let dy = self.y() - q.y();
        dx * dx + dy * dy
    }

    fn orient(&self, q: &Self, r: &Self) -> Orient {
        // Returns a **negative** value if ```self```, ```q``` and ```r``` occur in counterclockwise order (```r``` is to the left of the directed line ```self``` --> ```q```)
        // Returns a **positive** value if they occur in clockwise order(```r``` is to the right of the directed line ```self``` --> ```q```)
        // Returns zero is they are collinear
        match robust::orient2d(
            into_robust_coord(self),
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

    fn circumdelta(&self, b: &Self, c: &Self) -> (Self::Number, Self::Number) {
        let dx = b.x() - self.x();
        let dy = b.y() - self.y();
        let ex = c.x() - self.x();
        let ey = c.y() - self.y();

        let bl = dx * dx + dy * dy;
        let cl = ex * ex + ey * ey;
        let d = Self::Number::ONE_HALF / (dx * ey - dy * ex);

        let x = (ey * bl - dy * cl) * d;
        let y = (dx * cl - ex * bl) * d;
        (x, y)
    }

    fn circumradius2(&self, b: &Self, c: &Self) -> Self::Number {
        let (x, y) = self.circumdelta(b, c);
        x * x + y * y
    }

    fn circumcenter(&self, b: &Self, c: &Self) -> Self {
        let (x, y) = self.circumdelta(b, c);
        Self::new_point(self.x() + x, self.y() + y)
    }

    fn in_circle(&self, b: &Self, c: &Self, p: &Self) -> bool {
        let dx = self.x() - p.x();
        let dy = self.y() - p.y();
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

    fn nearly_equals(&self, p: &Self) -> bool {
        (self.x() - p.x()).abs() <= Self::Number::EPSILON
            && (self.y() - p.y()).abs() <= Self::Number::EPSILON
    }
}
