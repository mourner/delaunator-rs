#[macro_use]
extern crate criterion;

extern crate delaunator;
extern crate rand;

use criterion::{AxisScale, Criterion, ParameterizedBenchmark, PlotConfiguration, Throughput};
use delaunator::{triangulate, Point};
use rand::{Rng, SeedableRng, XorShiftRng};
use std::iter::repeat_with;

const COUNTS: &[usize] = &[100, 1000, 10_000, 100_000];

fn bench(c: &mut Criterion) {
    let mut rng = XorShiftRng::from_seed([0; 16]);

    let all_points: Vec<_> = repeat_with(|| rng.gen())
        .map(|(x, y)| Point { x, y })
        .take(*COUNTS.last().unwrap())
        .collect();

    c.bench(
        "triangulate",
        ParameterizedBenchmark::new(
            "triangulate",
            move |b, &&count| {
                let points = &all_points[..count];
                b.iter(move || triangulate(points))
            },
            COUNTS,
        )
        // need to override to a small sample size,
        // otherwise 100 000 elems takes 1-2 mins
        .sample_size(20)
            .plot_config(PlotConfiguration::default().summary_scale(AxisScale::Logarithmic))
            .throughput(|&&count| Throughput::Elements(count as _)),
    );
}

criterion_group!(benches, bench);
criterion_main!(benches);
