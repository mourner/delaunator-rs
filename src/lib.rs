use std::f64;
use std::fmt;

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
    fn dist2(&self, p: &Point) -> f64 {
        let dx = self.x - p.x;
        let dy = self.y - p.y;
        dx * dx + dy * dy
    }

    fn orient(&self, q: &Point, r: &Point) -> bool {
        (q.y - self.y) * (r.x - q.x) - (q.x - self.x) * (r.y - q.y) < 0.0
    }

    fn circumdelta(&self, b: &Point, c: &Point) -> (f64, f64) {
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

    fn circumradius2(&self, b: &Point, c: &Point) -> f64 {
        let (x, y) = self.circumdelta(b, c);
        x * x + y * y
    }

    fn circumcenter(&self, b: &Point, c: &Point) -> Point {
        let (x, y) = self.circumdelta(b, c);
        Point {
            x: self.x + x,
            y: self.y + y,
        }
    }

    fn in_circle(&self, b: &Point, c: &Point, p: &Point) -> bool {
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

    fn nearly_equals(&self, p: &Point) -> bool {
        (self.x - p.x).abs() <= f64::EPSILON && (self.y - p.y).abs() <= f64::EPSILON
    }
}

const EMPTY: usize = usize::max_value();

pub struct Triangulation {
    triangles: Vec<usize>,
    halfedges: Vec<usize>,
}

impl Triangulation {
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
        let a0 = a - a % 3;
        let ar = a0 + (a + 2) % 3;

        if b == EMPTY {
            return ar;
        }

        let b0 = b - b % 3;
        let al = a0 + (a + 1) % 3;
        let bl = b0 + (b + 2) % 3;

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

            let br = b0 + (b + 1) % 3;

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
    fn new(n: usize, center: Point, i0: usize, i1: usize, i2: usize, points: &[Point]) -> Hull {
        let hash_len = (n as f64).sqrt() as usize;

        let mut hull = Hull {
            prev: vec![0; n],            // edge to prev edge
            next: vec![0; n],            // edge to next edge
            tri: vec![0; n],             // edge to adjacent halfedge
            hash: vec![EMPTY; hash_len], // angular edge hash
            start: i0,
            center: center,
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

    fn find_visible_edge(&self, p: &Point, points: &[Point]) -> usize {
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
                return EMPTY;
            }
        }
        e
    }
}

pub fn triangulate(points: &[Point]) -> Triangulation {
    let n = points.len();
    let max_triangles = 2 * n - 5;

    // arrays that will store the triangulation graph
    let mut triangulation = Triangulation {
        triangles: Vec::with_capacity(max_triangles * 3),
        halfedges: Vec::with_capacity(max_triangles * 3),
    };

    // populate an array of point indices
    let mut ids: Vec<usize> = (0..n).collect();

    // calculate input data bbox center
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
    let bbox_center = Point {
        x: (min_x + max_x) / 2.0,
        y: (min_y + max_y) / 2.0,
    };

    // pick a seed point close to the center
    let mut min_dist = f64::INFINITY;
    let mut i0: usize = 0;
    for (i, p) in points.iter().enumerate() {
        let d = bbox_center.dist2(p);
        if d < min_dist {
            i0 = i;
            min_dist = d;
        }
    }
    let p0 = &points[i0];

    // find the point closest to the seed
    min_dist = f64::INFINITY;
    let mut i1: usize = 0;
    for (i, p) in points.iter().enumerate() {
        if i == i0 {
            continue;
        }
        let d = p0.dist2(p);
        if d < min_dist {
            i1 = i;
            min_dist = d;
        }
    }
    let mut p1 = &points[i1];

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
    let mut p2 = &points[i2];

    if min_radius == f64::INFINITY {
        panic!("No triangulation exists for this input");
    }

    // swap the order of the seed points for counter-clockwise orientation
    if p0.orient(p1, p2) {
        let i = i1;
        i1 = i2;
        i2 = i;
        let p = p1;
        p1 = p2;
        p2 = p;
    }

    // sort the points by distance from the seed triangle circumcenter
    let center = p0.circumcenter(p1, p2);

    ids.sort_unstable_by(|i, j| {
        let da = center.dist2(&points[*i]);
        let db = center.dist2(&points[*j]);
        da.partial_cmp(&db).unwrap()
    });

    let mut hull = Hull::new(n, center, i0, i1, i2, points);

    triangulation.add_triangle(i0, i1, i2, EMPTY, EMPTY, EMPTY);

    for (k, &i) in ids.iter().enumerate() {
        let p = &points[i];
        if k > 0 && p.nearly_equals(&points[ids[k - 1]]) {
            continue;
        }

        // skip seed triangle points
        if i == i0 || i == i1 || i == i2 {
            continue;
        }

        // find a visible edge on the convex hull using edge hash
        let mut e = hull.find_visible_edge(p, points);
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
            if p.orient(&points[n], &points[q]) {
                let t = triangulation.add_triangle(n, i, q, hull.tri[i], EMPTY, hull.tri[n]);
                hull.tri[i] = triangulation.legalize(t + 2, points, &mut hull);
                hull.next[n] = EMPTY; // mark as removed
                n = q;
            } else {
                break;
            }
        }

        // walk backward from the other side, adding more triangles and flipping
        loop {
            let q = hull.prev[e];
            if p.orient(&points[q], &points[e]) {
                let t = triangulation.add_triangle(q, i, e, EMPTY, hull.tri[e], hull.tri[q]);
                triangulation.legalize(t + 2, points, &mut hull);
                hull.tri[q] = t;
                hull.next[e] = EMPTY; // mark as removed
                e = q;
            } else {
                break;
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

    triangulation
}
