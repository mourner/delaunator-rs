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

    fn nearly_equals(&self, p: &Self) -> bool {
        (self.x - p.x).abs() <= EPSILON && (self.y - p.y).abs() <= EPSILON
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

    fn twin(&self, halfedge_id: usize) -> usize {
        self.halfedges[halfedge_id]
    }
    fn set_twin(&mut self, halfedge_id: usize, twin_id: usize) {
        if halfedge_id != EMPTY {
            self.halfedges[halfedge_id] = twin_id
        }
    }

    fn origin(&self, halfedge_id: usize) -> usize {
        self.triangles[halfedge_id]
    }
    fn set_origin(&mut self, halfedge_id: usize, point_id: usize) {
        self.triangles[halfedge_id] = point_id;
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

        self.set_twin(a, t);
        self.set_twin(b, t + 1);
        self.set_twin(c, t + 2);
        t
    }

    fn legalize(&mut self, a: usize, points: &[Point], hull: &mut Hull) -> usize {
        let b = self.twin(a);

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

        let p0 = self.origin(ar);
        let pr = self.origin(a);
        let pl = self.origin(al);
        let p1 = self.origin(bl);

        let illegal = (&points[p0]).in_circle(&points[pr], &points[pl], &points[p1]);
        if illegal {
            self.set_origin(a, p1);
            self.set_origin(b, p0);

            let hbl = self.twin(bl);
            let har = self.twin(ar);

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if hbl == EMPTY {
                hull.fix_halfedge(bl, a);
            }

            self.set_twin(a, hbl);
            self.set_twin(b, har);
            self.set_twin(ar, bl);

            self.set_twin(hbl, a);
            self.set_twin(har, b);
            self.set_twin(bl, ar);

            let br = next_halfedge(b);

            self.legalize(a, points, hull);
            return self.legalize(br, points, hull);
        }
        ar
    }
}

/// data structure for tracking the edges of the advancing convex hull
struct Hull {
    /// maps edge id to prev edge id
    prev: Vec<usize>,

    /// maps edge id to next edge id
    next: Vec<usize>,

    /// maps point id to outgoing halfedge id
    out: Vec<usize>,

    /// angular hull edge hash
    hash: Vec<usize>,

    /// starting point of the hull
    start: usize,

    /// center of the angular hash
    center: Point,
}

impl Hull {
    fn new(n: usize, center: Point, i0: usize, i1: usize, i2: usize, points: &[Point]) -> Self {
        let hash_len = (n as f64).sqrt() as usize;

        let mut hull = Self {
            prev: vec![0; n],
            next: vec![0; n],
            out: vec![0; n],
            hash: vec![EMPTY; hash_len],
            start: i0,
            center,
        };

        hull.set_next(i0, i1);
        hull.set_prev(i2, i1);
        hull.set_next(i1, i2);
        hull.set_prev(i0, i2);
        hull.set_next(i2, i0);
        hull.set_prev(i1, i0);

        hull.set_out(i0, 0);
        hull.set_out(i1, 1);
        hull.set_out(i2, 2);

        hull.hash_edge(&points[i0], i0);
        hull.hash_edge(&points[i1], i1);
        hull.hash_edge(&points[i2], i2);

        hull
    }

    fn out(&self, point_id: usize) -> usize {
        self.out[point_id]
    }
    fn set_out(&mut self, point_id: usize, halfedge_id: usize) {
        self.out[point_id] = halfedge_id;
    }

    fn prev(&self, point_id: usize) -> usize {
        self.prev[point_id]
    }
    fn set_prev(&mut self, point_id: usize, prev_point_id: usize) {
        self.prev[point_id] = prev_point_id;
    }

    fn next(&self, point_id: usize) -> usize {
        self.next[point_id]
    }
    fn set_next(&mut self, point_id: usize, next_point_id: usize) {
        self.next[point_id] = next_point_id;
    }

    fn remove(&mut self, point_id: usize) {
        self.set_next(point_id, EMPTY); // mark as removed
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
            if start != EMPTY && self.next(start) != EMPTY {
                break;
            }
        }
        start = self.prev(start);
        let mut e = start;

        while !p.orient(&points[e], &points[self.next(e)]) {
            e = self.next(e);
            if e == start {
                return (EMPTY, false);
            }
        }
        (e, e == start)
    }

    fn fix_halfedge(&mut self, old_id: usize, new_id: usize) {
        let mut e = self.start;
        loop {
            if self.out(e) == old_id {
                self.set_out(e, new_id);
                break;
            }
            e = self.next(e);
            if e == self.start {
                break;
            }
        }
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

    for (k, &(i, _)) in dists.iter().enumerate() {
        let p = &points[i];

        // skip near-duplicates
        if k > 0 && p.nearly_equals(&points[dists[k - 1].0]) {
            continue;
        }
        // skip seed triangle points
        if i == i0 || i == i1 || i == i2 {
            continue;
        }

        // find a visible edge on the convex hull using edge hash
        let (mut e, walk_back) = hull.find_visible_edge(p, points);
        if e == EMPTY {
            continue; // likely a near-duplicate point; skip it
        }

        // add the first triangle from the point
        let t = triangulation.add_triangle(e, i, hull.next(e), EMPTY, EMPTY, hull.out(e));

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        let out = triangulation.legalize(t + 2, points, &mut hull);
        hull.set_out(i, out);
        hull.set_out(e, t); // keep track of boundary triangles on the hull

        // walk forward through the hull, adding more triangles and flipping recursively
        let mut n = hull.next(e);
        loop {
            let q = hull.next(n);
            if !p.orient(&points[n], &points[q]) {
                break;
            }
            let t = triangulation.add_triangle(n, i, q, hull.out(i), EMPTY, hull.out(n));
            let out = triangulation.legalize(t + 2, points, &mut hull);;
            hull.set_out(i, out);
            hull.remove(n);
            n = q;
        }

        // walk backward from the other side, adding more triangles and flipping
        if walk_back {
            loop {
                let q = hull.prev[e];
                if !p.orient(&points[q], &points[e]) {
                    break;
                }
                let t = triangulation.add_triangle(q, i, e, EMPTY, hull.out(e), hull.out(q));
                triangulation.legalize(t + 2, points, &mut hull);
                hull.set_out(q, t);
                hull.remove(e);
                e = q;
            }
        }

        // update the hull indices
        hull.set_prev(i, e);
        hull.set_next(i, n);
        hull.set_prev(n, i);
        hull.set_next(e, i);
        hull.start = e;

        // save the two new edges in the hash table
        hull.hash_edge(p, i);
        hull.hash_edge(&points[e], e);
    }

    // expose hull as a vector of point indices
    let mut e = hull.start;
    loop {
        triangulation.hull.push(e);
        e = hull.next(e);
        if e == hull.start {
            break;
        }
    }

    triangulation.triangles.shrink_to_fit();
    triangulation.halfedges.shrink_to_fit();

    Some(triangulation)
}
