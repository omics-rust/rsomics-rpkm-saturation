use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_rpkm_saturation(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-rpkm-saturation");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    let bed = manifest.join("tests/golden/small.bed12");
    c.bench_function("rsomics-rpkm-saturation golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["-i", bam.to_str().unwrap(), "-r", bed.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_rpkm_saturation);
criterion_main!(benches);
