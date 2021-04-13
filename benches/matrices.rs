use criterion::{criterion_group, criterion_main, Criterion};
use gotoh::{ GotohMatrices, GotohParameters };

fn bench_matrices(c: &mut Criterion) {
    let mut group = c.benchmark_group("matrices");
    let count = 10000;
    let s = vec![b'A'; count];
    group.sample_size(10);
    group.bench_function(count.to_string(), |b| b.iter(
        || GotohMatrices::new(
            GotohParameters::new(1, -1, -3, -1),
            s.len(), s.len(),
            |x, y| s[x] == s[y])
    ));
    group.finish();
}

criterion_group!(benches, bench_matrices);
criterion_main!(benches);
