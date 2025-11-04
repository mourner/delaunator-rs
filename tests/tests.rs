use delaunator::{
    point_ext::{Orient, PointExt},
    triangulate, Number, Point, Triangulation, EMPTY,
};
use serde::de::DeserializeOwned;
use std::fmt::{Debug, Display};
use std::fs::File;

macro_rules! test_fixture {
    ($fixture_name:ident) => {
        #[test]
        fn $fixture_name() {
            let path = format!("tests/fixtures/{}.json", stringify!($fixture_name));
            validate(&load_fixture::<[f32; 2]>(&path));
            validate(&load_fixture::<[f64; 2]>(&path));
        }
    };
}

#[test]
fn basic() {
    validate(&load_fixture::<[f32; 2]>("tests/fixtures/basic.json"));
    validate(&load_fixture::<[f64; 2]>("tests/fixtures/basic.json"));
}

test_fixture!(robust2);
test_fixture!(robust3);
test_fixture!(robust4);
test_fixture!(robust5);
test_fixture!(robust6);
test_fixture!(issue10);
test_fixture!(ukraine);
test_fixture!(grid);

// issues from JS repo
test_fixture!(issue11js);
test_fixture!(issue13js);
test_fixture!(issue24js);
test_fixture!(issue43js);
test_fixture!(issue44js);

#[test]
fn robustness() {
    let points_f32 = load_fixture::<[f32; 2]>("tests/fixtures/robust1.json");

    validate(&points_f32);
    // f32 does not have enough precision to handle this scale
    // validate(&(scale_points(points_f32.iter().copied(), 1e-9)));
    validate(&(scale_points(points_f32.iter().copied(), 1e-2)));
    validate(&(scale_points(points_f32.iter().copied(), 100.0)));
    validate(&(scale_points(points_f32.iter().copied(), 1e9)));

    let points_f64 = load_fixture::<[f64; 2]>("tests/fixtures/robust1.json");

    validate(&points_f64);
    validate(&(scale_points(points_f64.iter().copied(), 1e-9)));
    validate(&(scale_points(points_f64.iter().copied(), 1e-2)));
    validate(&(scale_points(points_f64.iter().copied(), 100.0)));
    validate(&(scale_points(points_f64.iter().copied(), 1e9)));
}

#[test]
fn bad_input() {
    let mut points = vec![];
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (0 point)");
    assert!(halfedges.is_empty(), "Expected no edges (0 point)");
    assert!(hull.is_empty(), "Expected no hull (0 point)");

    points.push([0., 0.]);
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (1 point)");
    assert!(halfedges.is_empty(), "Expected no edges (1 point)");
    assert!(hull.len() == 1, "Expected single point on hull (1 point)");

    points.push([1., 0.]);
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (2 points)");
    assert!(halfedges.is_empty(), "Expected no edges (2 points)");
    assert!(hull.len() == 2, "Expected two points on hull (2 point)");
    assert!(
        hull.iter().enumerate().all(|(i, v)| i == *v),
        "Expected ordered hull points (2 point)"
    );

    points.push([2., 0.]);
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(
        triangles.is_empty(),
        "Expected no triangles (3 collinear points)"
    );
    assert!(
        halfedges.is_empty(),
        "Expected no edges (3 collinear points)"
    );
    assert_eq!(
        hull.len(),
        3,
        "Expected three points on hull (3 collinear points)"
    );
    assert!(
        hull.iter().enumerate().all(|(i, v)| i == *v),
        "Expected ordered hull points (3 collinear points)"
    );

    points.push([1., 1.]);
    validate(&points);
}

#[test]
fn unordered_collinear_points_input() {
    let points: Vec<_> = [10, 2, 4, 4, 1, 0, 3, 6, 8, 5, 7, 9]
        .iter()
        .map(|y| (0.0, *y as f64))
        .collect();
    let duplicated = 1;

    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(
        triangles.is_empty(),
        "Expected no triangles (unordered collinear points)"
    );
    assert!(
        halfedges.is_empty(),
        "Expected no edges (unordered collinear points)"
    );
    assert_eq!(
        hull.len(),
        points.len() - duplicated,
        "Expected all non-coincident points on hull (unordered collinear points)"
    );
    assert!(
        hull.iter()
            .enumerate()
            .all(|(i, v)| points[*v].y() == (i as f64)),
        "Expected ordered hull points (unordered collinear points)"
    );
}

#[test]
fn hull_collinear_issue24() {
    validate(&load_fixture::<[f32; 2]>("tests/fixtures/issue24.json"));

    let points = load_fixture::<[f64; 2]>("tests/fixtures/issue24.json");
    validate(&points);

    let t = triangulate(&points);
    assert_eq!(t.hull, &[0, 3, 2, 1], "Invalid hull");
}

#[test]
/// The test ensures that even when an invalid sequence of points is passed, there is no panic.
/// In this test, the output does not matter as long as an output is returned.
fn invalid_nan_sequence() {
    let points = vec![
        (-3.5, -1.5),
        (f64::NAN, f64::NAN),
        (f64::NAN, f64::NAN),
        (-3.5, -1.5),
    ];
    triangulate(&points);
}

fn scale_points<P: Point>(points: impl IntoIterator<Item = P>, scale: P::Number) -> Vec<P> {
    points
        .into_iter()
        .map(move |p| P::new_point(p.x() * scale, p.y() * scale))
        .collect()
}

fn load_fixture<P: Point>(path: &str) -> Vec<P>
where
    P::Number: DeserializeOwned,
{
    let file = File::open(path).unwrap();
    let u: Vec<(P::Number, P::Number)> = serde_json::from_reader(file).unwrap();
    u.iter().map(|p| Point::new_point(p.0, p.1)).collect()
}

/// make sure hull is convex and counter-clockwise (p1 is to the right of the directed line p0 --> p2)
/// in case of collinear points, make sure they are ordered (p1 between p0 and p2)
//  p-1                           p3
//   \                           ^
//    > p0 ---------------> p2 /
//              p1
fn assert_convex<N: Number + Into<f64>, P: Point<Number = N> + Debug>(p0: &P, p1: &P, p2: &P) {
    let l = p0.orient(p2, p1);
    assert!(matches!(l, Orient::Clockwise | Orient::Collinear), "p1 ({:?}) is to the left of the directed line p0 ({:?}) --> p2 ({:?}). Hull is not convex.", p1, p0, p2);

    if l == Orient::Collinear {
        // if p0, p1 and p2 are collinear, they must be ordered
        // that means that p1 - p0 = c * (p2 - p0), where c is (0..1) but not inclusive (linear combination)
        let c = ((p1.x() - p0.x()) / (p2.x() - p0.x())).max((p1.y() - p0.y()) / (p2.y() - p0.y()));
        assert!(c > N::ZERO, "incorrect ordering, found p1, p0, p2, expected p0 ({:?}), p1 ({:?}), p2 ({:?}). Invalid hull.", p0, p1, p2);
        assert!(c < N::ONE, "incorrect ordering, found p0, p2, p1, expected p0 ({:?}), p1 ({:?}), p2 ({:?}). Invalid hull.", p0, p1, p2);
    }
}

fn validate<N: Number + Display + Into<f64>>(points: &[impl Point<Number = N> + Debug]) {
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(points);

    // validate halfedges
    for (i, &h) in halfedges.iter().enumerate() {
        if h != EMPTY && halfedges[h] != i {
            panic!("Invalid halfedge connection");
        }
    }

    // validate triangulation
    let hull_area = {
        let mut hull_areas = Vec::new();

        for i in 0..hull.len() {
            let p0 = &points[hull[i]];
            let p1 = &points[hull[(i + 1) % hull.len()]];
            let p2 = &points[hull[(i + 2) % hull.len()]];
            assert_convex(p0, p1, p2);
            hull_areas.push((p1.x() - p0.x()) * (p1.y() + p0.y()));
        }

        sum(hull_areas)
    };
    let triangles_area = {
        let mut triangle_areas = Vec::new();
        let mut i = 0;
        while i < triangles.len() {
            let a = &points[triangles[i]];
            let b = &points[triangles[i + 1]];
            let c = &points[triangles[i + 2]];
            triangle_areas.push(
                ((b.y() - a.y()) * (c.x() - b.x()) - (b.x() - a.x()) * (c.y() - b.y())).abs(),
            );
            i += 3;
        }
        sum(triangle_areas)
    };

    let err = ((hull_area - triangles_area) / hull_area).abs();
    if err > N::EPSILON {
        panic!("Triangulation is broken: {} error", err);
    }
}

// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
fn sum<N: Number>(x: impl IntoIterator<Item = N>) -> N {
    let mut iter = x.into_iter();
    let mut sum = iter.next().unwrap_or(N::ZERO);
    let mut err = N::ZERO;
    for k in iter {
        let m = sum + k;
        err = err
            + if sum.abs() >= k.abs() {
                sum - m + k
            } else {
                k - m + sum
            };
        sum = m;
    }
    sum + err
}
