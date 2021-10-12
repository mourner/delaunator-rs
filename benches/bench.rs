#[macro_use]
extern crate criterion;

extern crate delaunator;
extern crate rand;

use criterion::{AxisScale, Criterion, ParameterizedBenchmark, PlotConfiguration};
use delaunator::{triangulate, Point};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::iter::repeat_with;

const COUNTS: &[usize] = &[100, 1000, 10_000, 100_000];

fn bench(c: &mut Criterion) {
    let mut rng: StdRng = StdRng::seed_from_u64(123);

    let all_points: Vec<_> = repeat_with(|| rng.gen())
        .map(|(x, y)| Point { x, y })
        .take(*COUNTS.last().unwrap())
        .collect();

    let bench = ParameterizedBenchmark::new(
        "triangulate",
        move |b, &&count| {
            let points = &all_points[..count];
            b.iter(move || triangulate(points))
        },
        COUNTS,
    );

    c.bench(
        "triangulate",
        bench
            .sample_size(20) // override to a small sample size, otherwise it takes too long
            .plot_config(PlotConfiguration::default().summary_scale(AxisScale::Logarithmic)),
    );
}

criterion_group!(benches, bench);
criterion_main!(benches);
