use delaunator::{triangulate, Point, Triangulation, EMPTY, EPSILON};
use std::f64;
use std::fs::File;

macro_rules! test_fixture {
    ($fixture_name:ident) => {
        #[test]
        fn $fixture_name() {
            let path = format!("tests/fixtures/{}.json", stringify!($fixture_name));
            let points = load_fixture(&path);
            validate(&points);
        }
    };
}

#[test]
fn basic() {
    validate(&load_fixture("tests/fixtures/basic.json"));
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

#[test]
fn robustness() {
    let points = load_fixture("tests/fixtures/robust1.json");

    validate(&points);
    validate(&(scale_points(&points, 1e-9)));
    validate(&(scale_points(&points, 1e-2)));
    validate(&(scale_points(&points, 100.0)));
    validate(&(scale_points(&points, 1e9)));
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

    points.push(Point { x: 0., y: 0. });
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (1 point)");
    assert!(halfedges.is_empty(), "Expected no edges (1 point)");
    assert!(hull.len() == 1, "Expected single point on hull (1 point)");

    points.push(Point { x: 1., y: 0. });
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (2 points)");
    assert!(halfedges.is_empty(), "Expected no edges (2 points)");
    assert!(hull.len() == 2, "Expected two points on hull (2 point)");
    assert!(hull.iter().enumerate().all(|(i, v)| i == *v), "Expected ordered hull points (2 point)");

    points.push(Point { x: 2., y: 0. });
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (3 collinear points)");
    assert!(halfedges.is_empty(), "Expected no edges (3 collinear points)");
    assert!(hull.len() == 3, "Expected three points on hull (3 collinear points)");
    assert!(hull.iter().enumerate().all(|(i, v)| i == *v), "Expected ordered hull points (3 collinear points)");

    points.push(Point { x: 1., y: 1. });
    validate(&points);
}

#[test]
fn unordered_collinear_points_input() {
    let points: Vec<Point> = [10, 2, 4, 4, 1, 0, 3, 6, 8, 5, 7, 9]
        .iter()
        .map(|y| Point { x: 0.0, y: *y as f64 })
        .collect();
    let duplicated = 1;

    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    assert!(triangles.is_empty(), "Expected no triangles (unordered collinear points)");
    assert!(halfedges.is_empty(), "Expected no edges (unordered collinear points)");
    assert!(hull.len() == points.len() - duplicated, "Expected all non-coincident points on hull (unordered collinear points)");
    assert!(hull.iter().enumerate().all(|(i, v)| points[*v].y == (i as f64)), "Expected ordered hull points (unordered collinear points)");
}

fn scale_points(points: &[Point], scale: f64) -> Vec<Point> {
    let scaled: Vec<Point> = points
        .iter()
        .map(|p| Point {
            x: p.x * scale,
            y: p.y * scale,
        }).collect();
    scaled
}

fn load_fixture(path: &str) -> Vec<Point> {
    let file = File::open(path).unwrap();
    let u: Vec<(f64, f64)> = serde_json::from_reader(file).unwrap();
    u.iter().map(|p| Point { x: p.0, y: p.1 }).collect()
}

fn validate(points: &[Point]) {
    let Triangulation {
        triangles,
        halfedges,
        hull,
    } = triangulate(&points);

    // validate halfedges
    for (i, &h) in halfedges.iter().enumerate() {
        if h != EMPTY && halfedges[h] != i {
            panic!("Invalid halfedge connection");
        }
    }

    // Currently failing for several fixtures: https://github.com/mourner/delaunator-rs/issues/18
    // validate connectivity
    // let mut seen = vec![false; points.len()];
    // for &i in triangles.iter() {
    //     seen[i] = true;
    //     seen[triangles[next_halfedge(i)]] = true;
    //     seen[triangles[next_halfedge(next_halfedge(i))]] = true;
    // }
    // let not_visited = seen.iter().enumerate().filter(|(_, visited)| !**visited).map(|(i,_)| i).collect::<Vec<usize>>();
    // assert_eq!(not_visited.len(), 0, "The following points are not part of the triangulation: {:#?}", not_visited);

    // validate triangulation
    let hull_area = {
        let mut hull_areas = Vec::new();
        let mut i = 0;
        let mut j = hull.len() - 1;
        while i < hull.len() {
            let p0 = &points[hull[j]];
            let p = &points[hull[i]];
            hull_areas.push((p.x - p0.x) * (p.y + p0.y));
            j = i;
            i += 1;
        }
        sum(&hull_areas)
    };
    let triangles_area = {
        let mut triangle_areas = Vec::new();
        let mut i = 0;
        while i < triangles.len() {
            let a = &points[triangles[i]];
            let b = &points[triangles[i + 1]];
            let c = &points[triangles[i + 2]];
            triangle_areas.push(((b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y)).abs());
            i += 3;
        }
        sum(&triangle_areas)
    };

    let err = ((hull_area - triangles_area) / hull_area).abs();
    if err > EPSILON {
        panic!("Triangulation is broken: {} error", err);
    }
}

// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
fn sum(x: &[f64]) -> f64 {
    let mut sum = x[0];
    let mut err: f64 = 0.0;
    for i in 1..x.len() {
        let k = x[i];
        let m = sum + k;
        err += if sum.abs() >= k.abs() {
            sum - m + k
        } else {
            k - m + sum
        };
        sum = m;
    }
    sum + err
}
