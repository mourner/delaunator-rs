# delaunator-rs

A very fast static 2D [Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) library for Rust.
A port of [Delaunator](https://github.com/mapbox/delaunator).

[![delaunator on Crates.io](https://meritbadge.herokuapp.com/delaunator?1)](https://crates.io/crates/delaunator)
[![Build Status](https://travis-ci.com/mourner/delaunator-rs.svg?branch=master)](https://travis-ci.com/mourner/delaunator-rs)

## [Documentation](https://docs.rs/delaunator)

## Example

```rust
use delaunator::{Point, triangulate};

let points = vec![
    Point { x: 0., y: 0. },
    Point { x: 1., y: 0. },
    Point { x: 1., y: 1. },
    Point { x: 0., y: 1. },
];

let result = triangulate(&points).expect("No triangulation exists.");

println!("{:?}", result.triangles); // [0, 2, 1, 0, 3, 2]
```

## Performance

Results for 3.1 GHz Intel Core i7 on a Macbook Pro 15'' (2017):

| points | time |
| ---: | ---: |
| 100 | 16.478µs |
| 1,000 | 277.64µs |
| 10,000 | 3.753ms |
| 100,000 | 63.627ms |
| 1,000,000 | 898.78ms |
| 10,000,000 | 11.857s |
