#[macro_use]
extern crate criterion;

extern crate delaunator;
extern crate rand;

use core::iter::repeat_with;
use criterion::{AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use delaunator::{triangulate, Point};
use rand::{rngs::StdRng, Rng, SeedableRng};

const COUNTS: &[usize] = &[100, 1000, 10_000, 100_000];

fn bench(c: &mut Criterion) {
    let mut rng: StdRng = StdRng::seed_from_u64(123);

    let all_points: Vec<_> = repeat_with(|| rng.random())
        .map(|(x, y)| Point { x, y })
        .take(*COUNTS.last().unwrap())
        .collect();

    let mut group = c.benchmark_group("triangulate");
    group.plot_config(PlotConfiguration::default().summary_scale(AxisScale::Logarithmic));

    for count in COUNTS {
        group.bench_with_input(BenchmarkId::from_parameter(count), count, |b, &count| {
            let points = &all_points[..count];
            b.iter(|| triangulate(points))
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
