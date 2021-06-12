extern crate delaunator;
extern crate serde_json;

use delaunator::{triangulate, Point, Triangulation, EMPTY, EPSILON};
use std::f64;
use std::fs::File;

#[test]
fn basic() {
    validate(&load_fixture("tests/fixtures/basic.json"));
}

#[test]
fn js_issues() {
    validate(&load_fixture("tests/fixtures/issue11.json"));
    validate(&load_fixture("tests/fixtures/issue13.json"));
    validate(&load_fixture("tests/fixtures/issue24.json"));
}

#[test]
fn robustness() {
    let points = load_fixture("tests/fixtures/robust1.json");

    validate(&points);
    validate(&(scale_points(&points, 1e-9)));
    validate(&(scale_points(&points, 1e-2)));
    validate(&(scale_points(&points, 100.0)));
    validate(&(scale_points(&points, 1e9)));

    validate(&load_fixture("tests/fixtures/robust2.json"));
    validate(&load_fixture("tests/fixtures/robust3.json"));
}

#[test]
fn bad_input() {
    let mut points = vec![Point { x: 0., y: 0. }];
    assert!(
        triangulate(&points).is_none(),
        "Expected empty triangulation (1 point)"
    );

    points.push(Point { x: 1., y: 0. });
    assert!(
        triangulate(&points).is_none(),
        "Expected empty triangulation (2 point)"
    );

    points.push(Point { x: 2., y: 0. });
    assert!(
        triangulate(&points).is_none(),
        "Expected empty triangulation (collinear points)"
    );

    points.push(Point { x: 1., y: 1. });
    validate(&points);
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
    } = triangulate(&points).expect("No triangulation exists for this input");

    // validate halfedges
    for (i, &h) in halfedges.iter().enumerate() {
        if h != EMPTY && halfedges[h] != i {
            panic!("Invalid halfedge connection");
        }
    }

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
