extern crate delaunator;
extern crate rand;

use std::iter::repeat_with;

const N: usize = 1_000_000;

fn main() {
    let points: Vec<_> = repeat_with(rand::random)
        .map(|(x, y)| delaunator::Point { x, y })
        .take(N)
        .collect();

    let now = std::time::Instant::now();
    let result = delaunator::triangulate(&points);
    let elapsed = now.elapsed();

    println!(
        "Triangulated {} points in {}.{}s.\nGenerated {} triangles. Convex hull size: {}",
        N,
        elapsed.as_secs(),
        elapsed.subsec_millis(),
        result.len(),
        result.hull.len()
    );
}
