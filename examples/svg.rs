use std::{env, fs::File, io::Write};
use delaunator::{EMPTY, Point, Triangulation, next_halfedge};
const CANVAS_SIZE: f64 = 800.;
const POINT_SIZE: usize = 4;
const LINE_WIDTH: usize = 1;
const HULL_COLOR: &str = "green";
const LINE_COLOR: &str = "blue";
const POINT_COLOR: &str = "black";
const HULL_POINT_COLOR: &str = "red";

/// Takes the first argument and use as path to load points data. If no argument provided, loads one of the test fixtures data file
/// Example: cargo run --example svg -- tests/fixtures/issue24.json
fn main() -> std::io::Result<()> {
    // load points from file
    let default_path = "tests/fixtures/robust4.json".to_string();
    let args = env::args().collect::<Vec<String>>();
    let path = args.get(1).unwrap_or(&default_path);
    let points: Vec<Point> = serde_json::from_reader::<_, Vec<(f64, f64)>>(File::open(path)?)?
        .iter().map(|p| Point { x: p.0, y: p.1 }).collect();

    // triangulate and scale points for display
    let triangulation = delaunator::triangulate(&points);
    println!("{:#?}", triangulation);
    let points = center_and_scale(&points, &triangulation);

    // generate SVG
    let contents = format!(
        r#"
<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">
<rect width="100%" height="100%" fill="white" />
    {circles}
    {lines}
</svg>"#,
        width = CANVAS_SIZE,
        height = CANVAS_SIZE,
        circles = render_point(&points, &triangulation),
        lines = (0..triangulation.triangles.len()).fold(String::new(), |acc, e| {
            if e > triangulation.halfedges[e] || triangulation.halfedges[e] == EMPTY {
                let start = &points[triangulation.triangles[e]];
                let end = &points[triangulation.triangles[next_halfedge(e)]];
                let color = if triangulation.halfedges[e] == EMPTY { HULL_COLOR } else { LINE_COLOR };
                acc + &format!(r#"<line x1="{x0}" y1="{y0}" x2="{x1}" y2="{y1}" style="stroke:{color};stroke-width:{width}" />"#, x0 = start.x, y0 = start.y, x1=end.x, y1=end.y, width = LINE_WIDTH, color = color)
            } else {
                acc
            }
        })
    );
    File::create("triangulation.svg")?
        .write_all(contents.as_bytes())
}


/// Finds the center point and farthest point from it, then generates a new vector of
/// scaled and offset points such that they fit between [0..SIZE]
fn center_and_scale(points: &Vec<Point>, t: &Triangulation) -> Vec<Point> {
    let center = &points[*t.triangles.get(0).unwrap_or(&0)];
    let farthest_distance = points.iter().map(|p| {
        let (x, y) = (center.x - p.x, center.y - p.y);
        x*x + y*y
    }).reduce(f64::max).unwrap().sqrt();
    let scale = CANVAS_SIZE / (farthest_distance * 2.0);
    let offset = ((CANVAS_SIZE / 2.0) - (scale * center.x), (CANVAS_SIZE / 2.0) - (scale * center.y));
    points.iter().map(|p| Point { x: scale * p.x + offset.0, y: scale * p.y + offset.1 }).collect()
}


fn render_point(points: &[Point], triangulation: &Triangulation) -> String {
    let mut circles = points.iter().enumerate().fold(String::new(), |acc, (i, p)| {
        let color = if triangulation.hull.contains(&i) { HULL_POINT_COLOR } else { POINT_COLOR };
        acc + &format!(r#"<circle cx="{x}" cy="{y}" r="{size}" fill="{color}"/>"#, x = p.x, y = p.y, size = POINT_SIZE, color = color)
    });

    // show ids for points if input is relatively small
    if points.len() < 100 {
        circles = points.iter().enumerate().fold(circles, |acc, (i, p)| {
            acc + &format!(r#"<text x="{x}" y="{y}" font-size="20" fill="black">{i}</text>"#, i = i, x = p.x + 10., y = p.y - 5.)
        })
    }

    circles
}