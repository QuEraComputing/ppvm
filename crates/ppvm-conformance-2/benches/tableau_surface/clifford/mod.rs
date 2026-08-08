// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;

use super::*;

mod batches;
mod blocks;
mod gates;

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/clifford");
    let (old_b, new_b) = prepared_bare(96);
    let (old_g, new_g) = prepared_gen(96);

    gates::bench(&mut group, &old_b, &new_b, &old_g, &new_g);
    batches::bench(&mut group, &old_b, &new_b, &old_g, &new_g);
    blocks::bench(&mut group, &old_b, &new_b, &old_g, &new_g);
    blocks::width_sweep(&mut group);
    group.finish();
}
