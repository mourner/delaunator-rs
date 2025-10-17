use core::iter::repeat_with;
use criterion::{criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion, PlotConfiguration};
use itertools::Itertools;
use delaunator::{triangulate, Point};
use rand::{rngs::StdRng, Rng, SeedableRng};

const COUNTS: &[usize] = &[100, 1000, 10_000, 100_000];

#[derive(Copy, Clone, Debug)]
enum PointType {
    DelaunatorPoint,
    RobustCoord,
}

fn bench(c: &mut Criterion) {
    let mut rng: StdRng = StdRng::seed_from_u64(123);

    let all_points: Vec<_> = repeat_with(|| rng.random())
        .map(|(x, y)| Point { x, y })
        .take(*COUNTS.last().unwrap())
        .collect();
    let other_point_type = all_points.iter().copied().map(robust::Coord::from).collect::<Vec<_>>();

    let mut group = c.benchmark_group("triangulate");
    group.plot_config(PlotConfiguration::default().summary_scale(AxisScale::Logarithmic));

    for (count, point_type) in COUNTS.iter().cartesian_product([PointType::DelaunatorPoint, PointType::RobustCoord]) {
        group.bench_with_input(BenchmarkId::from_parameter(format!("{count}-{point_type:?}")), &(count, point_type), |b, (&count, point_type)| {
            match point_type {
                PointType::DelaunatorPoint => {
                    let points = &all_points[..count];
                    b.iter(|| triangulate(points))
                }
                PointType::RobustCoord => {
                    let points = &other_point_type[..count];
                    b.iter(|| triangulate(points))
                }
            }
        });
    }
}

criterion_group!(benches, bench);
criterion_main!(benches);
