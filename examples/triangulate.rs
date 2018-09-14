extern crate delaunator;
extern crate rand;

fn main() {
    let n = 1000_000;
    let mut points = Vec::with_capacity(n);
    for _i in 0..n {
        points.push(delaunator::Point {
            x: rand::random(),
            y: rand::random(),
        });
    }

    let now = std::time::Instant::now();
    let result = delaunator::triangulate(&points);
    let elapsed = now.elapsed();

    println!(
        "Triangulated {} points in {}.{}s.\nGenerated {} triangles. Convex hull size: {}",
        n,
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        result.len(),
        result.hull.len()
    );
}
