/*!
A very fast 2D [Delaunay Triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) library for Rust.
A port of [Delaunator](https://github.com/mapbox/delaunator).

# Example

```rust
use delaunator::{Point, triangulate};

let points = vec![
    Point { x: 0., y: 0. },
    Point { x: 1., y: 0. },
    Point { x: 1., y: 1. },
    Point { x: 0., y: 1. },
];

let result = triangulate(&points).expect("No triangulation exists.");
println!("{:?}", result.triangles); // [0, 2, 1, 0, 3, 2]
```
*/

use std::{f64, fmt};

/// Near-duplicate points (where both `x` and `y` only differ within this value)
/// will not be included in the triangulation for robustness.
pub const EPSILON: f64 = f64::EPSILON * 2.0;

/// Represents a 2D point in the input vector.
#[derive(Clone, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl fmt::Debug for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{}, {}]", self.x, self.y)
    }
}

impl Point {
    fn dist2(&self, p: &Self) -> f64 {
        let dx = self.x - p.x;
        let dy = self.y - p.y;
        dx * dx + dy * dy
    }

    fn orient(&self, q: &Self, r: &Self) -> bool {
        (q.y - self.y) * (r.x - q.x) - (q.x - self.x) * (r.y - q.y) < 0.0
    }

    fn circumdelta(&self, b: &Self, c: &Self) -> (f64, f64) {
        let dx = b.x - self.x;
        let dy = b.y - self.y;
        let ex = c.x - self.x;
        let ey = c.y - self.y;

        let bl = dx * dx + dy * dy;
        let cl = ex * ex + ey * ey;
        let d = 0.5 / (dx * ey - dy * ex);

        let x = (ey * bl - dy * cl) * d;
        let y = (dx * cl - ex * bl) * d;
        (x, y)
    }

    fn circumradius2(&self, b: &Self, c: &Self) -> f64 {
        let (x, y) = self.circumdelta(b, c);
        x * x + y * y
    }

    fn circumcenter(&self, b: &Self, c: &Self) -> Self {
        let (x, y) = self.circumdelta(b, c);
        Self {
            x: self.x + x,
            y: self.y + y,
        }
    }

    fn in_circle(&self, b: &Self, c: &Self, p: &Self) -> bool {
        let dx = self.x - p.x;
        let dy = self.y - p.y;
        let ex = b.x - p.x;
        let ey = b.y - p.y;
        let fx = c.x - p.x;
        let fy = c.y - p.y;

        let ap = dx * dx + dy * dy;
        let bp = ex * ex + ey * ey;
        let cp = fx * fx + fy * fy;

        dx * (ey * cp - bp * fy) - dy * (ex * cp - bp * fx) + ap * (ex * fy - ey * fx) < 0.0
    }
}

/// Represents the area outside of the triangulation.
/// Halfedges on the convex hull (which don't have an adjacent halfedge)
/// will have this value.
pub const EMPTY: usize = usize::max_value();

/// Next halfedge in a triangle.
pub fn next_halfedge(i: usize) -> usize {
    if i % 3 == 2 {
        i - 2
    } else {
        i + 1
    }
}

/// Previous halfedge in a triangle.
pub fn prev_halfedge(i: usize) -> usize {
    if i % 3 == 0 {
        i + 2
    } else {
        i - 1
    }
}

/// Result of the Delaunay triangulation.
pub struct Triangulation {
    /// A vector of point indices where each triple represents a Delaunay triangle.
    /// All triangles are directed counter-clockwise.
    pub triangles: Vec<usize>,

    /// A vector of adjacent halfedge indices that allows traversing the triangulation graph.
    ///
    /// `i`-th half-edge in the array corresponds to vertex `triangles[i]`
    /// the half-edge is coming from. `halfedges[i]` is the index of a twin half-edge
    /// in an adjacent triangle (or `EMPTY` for outer half-edges on the convex hull).
    pub halfedges: Vec<usize>,

    /// A vector of indices that reference points on the convex hull of the triangulation,
    /// counter-clockwise.
    pub hull: Vec<usize>,
}

impl Triangulation {
    fn new(n: usize) -> Self {
        let max_triangles = 2 * n - 5;
        Self {
            triangles: Vec::with_capacity(max_triangles * 3),
            halfedges: Vec::with_capacity(max_triangles * 3),
            hull: Vec::new(),
        }
    }

    /// The number of triangles in the triangulation.
    pub fn len(&self) -> usize {
        (self.triangles.len() / 3)
    }

    fn add_triangle(
        &mut self,
        i0: usize,
        i1: usize,
        i2: usize,
        a: usize,
        b: usize,
        c: usize,
    ) -> usize {
        let t = self.triangles.len();

        self.triangles.push(i0);
        self.triangles.push(i1);
        self.triangles.push(i2);

        self.halfedges.push(a);
        self.halfedges.push(b);
        self.halfedges.push(c);

        if a != EMPTY {
            self.halfedges[a] = t;
        }
        if b != EMPTY {
            self.halfedges[b] = t + 1;
        }
        if c != EMPTY {
            self.halfedges[c] = t + 2;
        }

        t
    }

    fn legalize(&mut self, a: usize, points: &[Point], hull: &mut Hull) -> usize {
        let b = self.halfedges[a];

        // if the pair of triangles doesn't satisfy the Delaunay condition
        // (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
        // then do the same check/flip recursively for the new pair of triangles
        //
        //           pl                    pl
        //          /||\                  /  \
        //       al/ || \bl            al/    \a
        //        /  ||  \              /      \
        //       /  a||b  \    flip    /___ar___\
        //     p0\   ||   /p1   =>   p0\---bl---/p1
        //        \  ||  /              \      /
        //       ar\ || /br             b\    /br
        //          \||/                  \  /
        //           pr                    pr
        //
        let ar = prev_halfedge(a);

        if b == EMPTY {
            return ar;
        }

        let al = next_halfedge(a);
        let bl = prev_halfedge(b);

        let p0 = self.triangles[ar];
        let pr = self.triangles[a];
        let pl = self.triangles[al];
        let p1 = self.triangles[bl];

        let illegal = (&points[p0]).in_circle(&points[pr], &points[pl], &points[p1]);
        if illegal {
            self.triangles[a] = p1;
            self.triangles[b] = p0;

            let hbl = self.halfedges[bl];
            let har = self.halfedges[ar];

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if hbl == EMPTY {
                let mut e = hull.start;
                loop {
                    if hull.tri[e] == bl {
                        hull.tri[e] = a;
                        break;
                    }
                    e = hull.next[e];
                    if e == hull.start {
                        break;
                    }
                }
            }

            self.halfedges[a] = hbl;
            self.halfedges[b] = har;
            self.halfedges[ar] = bl;

            if hbl != EMPTY {
                self.halfedges[hbl] = a;
            }
            if har != EMPTY {
                self.halfedges[har] = b;
            }
            if bl != EMPTY {
                self.halfedges[bl] = ar;
            }

            let br = next_halfedge(b);

            self.legalize(a, points, hull);
            return self.legalize(br, points, hull);
        }
        ar
    }
}

// data structure for tracking the edges of the advancing convex hull
struct Hull {
    prev: Vec<usize>,
    next: Vec<usize>,
    tri: Vec<usize>,
    hash: Vec<usize>,
    start: usize,
    center: Point,
}

impl Hull {
    fn new(n: usize, center: Point, i0: usize, i1: usize, i2: usize, points: &[Point]) -> Self {
        let hash_len = (n as f64).sqrt() as usize;

        let mut hull = Self {
            prev: vec![0; n],            // edge to prev edge
            next: vec![0; n],            // edge to next edge
            tri: vec![0; n],             // edge to adjacent halfedge
            hash: vec![EMPTY; hash_len], // angular edge hash
            start: i0,
            center,
        };

        hull.next[i0] = i1;
        hull.prev[i2] = i1;
        hull.next[i1] = i2;
        hull.prev[i0] = i2;
        hull.next[i2] = i0;
        hull.prev[i1] = i0;

        hull.tri[i0] = 0;
        hull.tri[i1] = 1;
        hull.tri[i2] = 2;

        hull.hash_edge(&points[i0], i0);
        hull.hash_edge(&points[i1], i1);
        hull.hash_edge(&points[i2], i2);

        hull
    }

    fn hash_key(&self, p: &Point) -> usize {
        let dx = p.x - self.center.x;
        let dy = p.y - self.center.y;

        let p = dx / (dx.abs() + dy.abs());
        let a = (if dy > 0.0 { 3.0 - p } else { 1.0 + p }) / 4.0; // [0..1]

        let len = self.hash.len();
        (((len as f64) * a).floor() as usize) % len
    }

    fn hash_edge(&mut self, p: &Point, i: usize) {
        let key = self.hash_key(p);
        self.hash[key] = i;
    }

    fn find_visible_edge(&self, p: &Point, points: &[Point]) -> (usize, bool) {
        let mut start: usize = 0;
        let key = self.hash_key(p);
        let len = self.hash.len();
        for j in 0..len {
            start = self.hash[(key + j) % len];
            if start != EMPTY && self.next[start] != EMPTY {
                break;
            }
        }
        start = self.prev[start];
        let mut e = start;

        while !p.orient(&points[e], &points[self.next[e]]) {
            e = self.next[e];
            if e == start {
                return (EMPTY, false);
            }
        }
        (e, e == start)
    }
}

fn calc_bbox_center(points: &[Point]) -> Point {
    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for p in points.iter() {
        min_x = min_x.min(p.x);
        min_y = min_y.min(p.y);
        max_x = max_x.max(p.x);
        max_y = max_y.max(p.y);
    }
    Point {
        x: (min_x + max_x) / 2.0,
        y: (min_y + max_y) / 2.0,
    }
}

fn find_closest_point(points: &[Point], p0: &Point) -> Option<usize> {
    let mut min_dist = f64::INFINITY;
    let mut k: usize = 0;
    for (i, p) in points.iter().enumerate() {
        let d = p0.dist2(p);
        if d > 0.0 && d < min_dist {
            k = i;
            min_dist = d;
        }
    }
    if min_dist == f64::INFINITY {
        None
    } else {
        Some(k)
    }
}

fn find_seed_triangle(points: &[Point]) -> Option<(usize, usize, usize)> {
    // pick a seed point close to the center
    let bbox_center = calc_bbox_center(points);
    let i0 = find_closest_point(points, &bbox_center)?;
    let p0 = &points[i0];

    // find the point closest to the seed
    let i1 = find_closest_point(points, p0)?;
    let p1 = &points[i1];

    // find the third point which forms the smallest circumcircle with the first two
    let mut min_radius = f64::INFINITY;
    let mut i2: usize = 0;
    for (i, p) in points.iter().enumerate() {
        if i == i0 || i == i1 {
            continue;
        }
        let r = p0.circumradius2(p1, p);
        if r < min_radius {
            i2 = i;
            min_radius = r;
        }
    }

    if min_radius == f64::INFINITY {
        None
    } else {
        // swap the order of the seed points for counter-clockwise orientation
        Some(if p0.orient(p1, &points[i2]) {
            (i0, i2, i1)
        } else {
            (i0, i1, i2)
        })
    }
}

/// Triangulate a set of 2D points.
/// Returns `None` if no triangulation exists for the input (e.g. all points are collinear).
pub fn triangulate(points: &[Point]) -> Option<Triangulation> {
    let n = points.len();

    let (i0, i1, i2) = find_seed_triangle(points)?;
    let center = (&points[i0]).circumcenter(&points[i1], &points[i2]);

    let mut triangulation = Triangulation::new(n);
    triangulation.add_triangle(i0, i1, i2, EMPTY, EMPTY, EMPTY);

    // sort the points by distance from the seed triangle circumcenter
    let mut dists: Vec<_> = points
        .iter()
        .enumerate()
        .map(|(i, point)| (i, center.dist2(point)))
        .collect();

    dists.sort_unstable_by(|&(_, da), &(_, db)| da.partial_cmp(&db).unwrap());

    let mut hull = Hull::new(n, center, i0, i1, i2, points);

    let mut prev_dist = f64::INFINITY;

    for (i, d) in dists {
        let p = &points[i];

        // skip seed triangle points
        if i == i0 || i == i1 || i == i2 {
            continue;
        }

        // skip near-duplicates
        if (d - prev_dist).abs() <= EPSILON {
            continue;
        }

        prev_dist = d;

        // find a visible edge on the convex hull using edge hash
        let (mut e, walk_back) = hull.find_visible_edge(p, points);
        if e == EMPTY {
            continue; // likely a near-duplicate point; skip it
        }

        // add the first triangle from the point
        let t = triangulation.add_triangle(e, i, hull.next[e], EMPTY, EMPTY, hull.tri[e]);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
        hull.tri[e] = t; // keep track of boundary triangles on the hull

        // walk forward through the hull, adding more triangles and flipping recursively
        let mut n = hull.next[e];
        loop {
            let q = hull.next[n];
            if !p.orient(&points[n], &points[q]) {
                break;
            }
            let t = triangulation.add_triangle(n, i, q, hull.tri[i], EMPTY, hull.tri[n]);
            hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
            hull.next[n] = EMPTY; // mark as removed
            n = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if walk_back {
            loop {
                let q = hull.prev[e];
                if !p.orient(&points[q], &points[e]) {
                    break;
                }
                let t = triangulation.add_triangle(q, i, e, EMPTY, hull.tri[e], hull.tri[q]);
                triangulation.legalize(t + 2, points, &mut hull);
                hull.tri[q] = t;
                hull.next[e] = EMPTY; // mark as removed
                e = q;
            }
        }

        // update the hull indices
        hull.prev[i] = e;
        hull.next[i] = n;
        hull.prev[n] = i;
        hull.next[e] = i;
        hull.start = e;

        // save the two new edges in the hash table
        hull.hash_edge(p, i);
        hull.hash_edge(&points[e], e);
    }

    // expose hull as a vector of point indices
    let mut e = hull.start;
    loop {
        triangulation.hull.push(e);
        e = hull.next[e];
        if e == hull.start {
            break;
        }
    }

    triangulation.triangles.shrink_to_fit();
    triangulation.halfedges.shrink_to_fit();

    Some(triangulation)
}
