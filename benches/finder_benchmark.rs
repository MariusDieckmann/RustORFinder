use criterion::{criterion_group, criterion_main, Criterion};
use rustyorffinder::finder::threaded_finder::ThreadedFinder;

fn criterion_benchmark(c: &mut Criterion) {
    let sequence = include_str!("../resources/sequences/NC_011604.1_normal.fasta");

    let threaded_finder = ThreadedFinder::new(sequence.to_string());

    c.bench_function("threaded_orf_finder", |b| {
        b.iter(|| threaded_finder.find_orfs())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
