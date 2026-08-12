// SPDX-FileCopyrightText: 2026 The PPVM Authors
// SPDX-License-Identifier: Apache-2.0

use criterion::Criterion;

use super::*;

macro_rules! channel {
    ($group:expr, $name:expr, $trait:ident, $op:ident, $old_b:expr, $new_b:expr,
     $old_g:expr, $new_g:expr $(, $arg:expr)* $(,)?) => {{
        let (mut ob, mut nb) = (($old_b).clone(), ($new_b).clone());
        ppvm_traits::traits::$trait::$op(&mut ob $(, $arg)*);
        ppvm_traits_2::$trait::$op(
            &mut nb
            $(, $arg)*,
            &mut ppvm_conformance_2::analytic_rng(),
        );
        assert_bare_eq(&ob, &nb);
        let (mut og, mut ng) = (($old_g).clone(), ($new_g).clone());
        ppvm_traits::traits::$trait::$op(&mut og $(, $arg)*);
        ppvm_traits_2::$trait::$op(
            &mut ng
            $(, $arg)*,
            &mut ppvm_conformance_2::analytic_rng(),
        );
        assert_gen_eq(&og, &ng);
        bench_mut_pair!($group, concat!("bare/", $name), $old_b, $new_b,
            |t: &mut OldBare| ppvm_traits::traits::$trait::$op(t $(, $arg)*),
            |t: &mut NewBare| ppvm_traits_2::$trait::$op(
                t
                $(, $arg)*,
                &mut ppvm_conformance_2::analytic_rng(),
            ));
        bench_mut_pair!($group, concat!("generalized/", $name), $old_g, $new_g,
            |t: &mut OldGen| ppvm_traits::traits::$trait::$op(t $(, $arg)*),
            |t: &mut NewGen| ppvm_traits_2::$trait::$op(
                t
                $(, $arg)*,
                &mut ppvm_conformance_2::analytic_rng(),
            ));
    }};
}

macro_rules! loss {
    ($group:expr, $name:expr, $trait:ident, $op:ident, $old:expr, $new:expr
     $(, $arg:expr)* $(,)?) => {{
        let (mut oc, mut nc) = (($old).clone(), ($new).clone());
        ppvm_traits::traits::$trait::$op(&mut oc $(, $arg)*);
        ppvm_traits_2::$trait::$op(
            &mut nc
            $(, $arg)*,
            &mut ppvm_conformance_2::analytic_rng(),
        );
        assert_gen_eq(&oc, &nc);
        bench_mut_pair!($group, concat!("generalized/", $name), $old, $new,
            |t: &mut OldGen| ppvm_traits::traits::$trait::$op(t $(, $arg)*),
            |t: &mut NewGen| ppvm_traits_2::$trait::$op(
                t
                $(, $arg)*,
                &mut ppvm_conformance_2::analytic_rng(),
            ));
    }};
}

pub fn bench(c: &mut Criterion) {
    let mut group = c.benchmark_group("tableau-surface/noise");
    let (old_b, new_b) = prepared_bare(96);
    let (old_g, new_g) = prepared_gen(96);
    let qs = [0usize, 2, 4, 66];
    let pairs = [(0usize, 1usize), (2, 3), (64, 65)];
    let pauli = [0.1f64, 0.2, 0.3];
    let two = [0.02f64; 15];

    channel!(
        group,
        "pauli_error",
        PauliError,
        pauli_error,
        old_b,
        new_b,
        old_g,
        new_g,
        0,
        pauli
    );
    channel!(
        group, "x_error", PauliError, x_error, old_b, new_b, old_g, new_g, 0, 0.3
    );
    channel!(
        group, "y_error", PauliError, y_error, old_b, new_b, old_g, new_g, 0, 0.3
    );
    channel!(
        group, "z_error", PauliError, z_error, old_b, new_b, old_g, new_g, 0, 0.3
    );
    channel!(
        group,
        "pauli_error_many",
        PauliError,
        pauli_error_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &qs,
        pauli
    );
    channel!(
        group,
        "x_error_many",
        PauliError,
        x_error_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &qs,
        0.3
    );
    channel!(
        group,
        "y_error_many",
        PauliError,
        y_error_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &qs,
        0.3
    );
    channel!(
        group,
        "z_error_many",
        PauliError,
        z_error_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &qs,
        0.3
    );
    channel!(
        group,
        "two_qubit_pauli_error",
        TwoQubitPauliError,
        two_qubit_pauli_error,
        old_b,
        new_b,
        old_g,
        new_g,
        0,
        65,
        two
    );
    channel!(
        group,
        "two_qubit_pauli_error_many",
        TwoQubitPauliError,
        two_qubit_pauli_error_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &pairs,
        two
    );
    channel!(
        group,
        "depolarize1",
        Depolarizing,
        depolarize1,
        old_b,
        new_b,
        old_g,
        new_g,
        0,
        0.3
    );
    channel!(
        group,
        "depolarize1_many",
        Depolarizing,
        depolarize1_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &qs,
        0.3
    );
    channel!(
        group,
        "depolarize2",
        Depolarizing2,
        depolarize2,
        old_b,
        new_b,
        old_g,
        new_g,
        0,
        65,
        0.3
    );
    channel!(
        group,
        "depolarize2_many",
        Depolarizing2,
        depolarize2_many,
        old_b,
        new_b,
        old_g,
        new_g,
        &pairs,
        0.3
    );

    loss!(
        group,
        "loss_channel",
        LossChannel,
        loss_channel,
        old_g,
        new_g,
        0,
        0.3
    );
    loss!(
        group,
        "asymmetric_loss_channel",
        AsymmetricLossChannel,
        asymmetric_loss_channel,
        old_g,
        new_g,
        0,
        0.2,
        0.4
    );
    loss!(
        group,
        "correlated_loss_channel",
        CorrelatedLossChannel,
        correlated_loss_channel,
        old_g,
        new_g,
        0,
        65,
        [0.1, 0.2, 0.3]
    );

    let (mut old_lost, mut new_lost) = (old_g.clone(), new_g.clone());
    ppvm_traits::traits::LossChannel::loss_channel(&mut old_lost, 0, 1.0);
    ppvm_traits_2::LossChannel::loss_channel(
        &mut new_lost,
        0,
        1.0,
        &mut ppvm_conformance_2::analytic_rng(),
    );
    let (mut old_check, mut new_check) = (old_lost.clone(), new_lost.clone());
    ppvm_traits::traits::ResetLossChannel::reset_loss_channel(&mut old_check, 0);
    ppvm_traits_2::ResetLossChannel::reset_loss_channel(&mut new_check, 0);
    assert_gen_eq(&old_check, &new_check);
    bench_mut_pair!(
        group,
        "generalized/reset_loss_channel",
        old_lost,
        new_lost,
        |t: &mut OldGen| ppvm_traits::traits::ResetLossChannel::reset_loss_channel(t, 0),
        |t: &mut NewGen| ppvm_traits_2::ResetLossChannel::reset_loss_channel(t, 0)
    );
    group.finish();
}
