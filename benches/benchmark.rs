use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fusion_neutron_utils::{scale_temperature_units_to_kev, tt_spectrum_per_ev};

fn benchmark_temperature_conversion(c: &mut Criterion) {
    c.bench_function("scale_temperature_units_to_kev", |b| {
        b.iter(|| {
            let _ = scale_temperature_units_to_kev(
                black_box(30_000.0),
                black_box(Some("eV")),
            );
        });
    });
}

fn benchmark_tt_spectrum(c: &mut Criterion) {
    let energies: Vec<f64> = (0..100_000).map(|i| i as f64 / 10.0).collect();
    c.bench_function("tt_spectrum_per_ev_100k", |b| {
        b.iter(|| {
            let _ = tt_spectrum_per_ev(black_box(&energies), black_box(10_000.0));
        });
    });
}

criterion_group!(benches, benchmark_temperature_conversion, benchmark_tt_spectrum);
criterion_main!(benches);
