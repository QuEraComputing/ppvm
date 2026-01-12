#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline(always)]
pub(crate) fn xor_phase_with(phase: &mut [u64], bits: &[u64]) {
    debug_assert_eq!(phase.len(), bits.len());
    if phase.len() < 4 {
        xor_phase_with_scalar(phase, bits);
    } else {
        unsafe { xor_with_avx2(phase, bits) };
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub(crate) fn xor_phase_with(phase: &mut [u64], bits: &[u64]) {
    debug_assert_eq!(phase.len(), bits.len());
    if phase.len() < 2 {
        xor_phase_with_scalar(phase, bits);
    } else {
        unsafe { xor_with_neon(phase, bits) };
    }
}

#[cfg(all(
    target_arch = "x86_64",
    not(target_feature = "avx2"),
    target_feature = "sse2"
))]
#[inline(always)]
pub(crate) fn xor_phase_with(phase: &mut [u64], bits: &[u64]) {
    debug_assert_eq!(phase.len(), bits.len());
    if phase.len() < 2 {
        xor_phase_with_scalar(phase, bits);
    } else {
        unsafe { xor_with_sse2(phase, bits) };
    }
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon"),
    all(
        target_arch = "x86_64",
        not(target_feature = "avx2"),
        target_feature = "sse2"
    )
)))]
#[inline(always)]
pub(crate) fn xor_phase_with(phase: &mut [u64], bits: &[u64]) {
    debug_assert_eq!(phase.len(), bits.len());
    xor_phase_with_scalar(phase, bits);
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
#[inline(always)]
pub(crate) fn xor_phase_with_xor(phase: &mut [u64], left: &[u64], right: &[u64]) {
    debug_assert_eq!(phase.len(), left.len());
    debug_assert_eq!(left.len(), right.len());
    if phase.len() < 4 {
        xor_phase_with_xor_scalar(phase, left, right);
    } else {
        unsafe { xor_with_xor_avx2(phase, left, right) };
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
pub(crate) fn xor_phase_with_xor(phase: &mut [u64], left: &[u64], right: &[u64]) {
    debug_assert_eq!(phase.len(), left.len());
    debug_assert_eq!(left.len(), right.len());
    if phase.len() < 2 {
        xor_phase_with_xor_scalar(phase, left, right);
    } else {
        unsafe { xor_with_xor_neon(phase, left, right) };
    }
}

#[cfg(all(
    target_arch = "x86_64",
    not(target_feature = "avx2"),
    target_feature = "sse2"
))]
#[inline(always)]
pub(crate) fn xor_phase_with_xor(phase: &mut [u64], left: &[u64], right: &[u64]) {
    debug_assert_eq!(phase.len(), left.len());
    debug_assert_eq!(left.len(), right.len());
    if phase.len() < 2 {
        xor_phase_with_xor_scalar(phase, left, right);
    } else {
        unsafe { xor_with_xor_sse2(phase, left, right) };
    }
}

#[cfg(not(any(
    all(target_arch = "x86_64", target_feature = "avx2"),
    all(target_arch = "aarch64", target_feature = "neon"),
    all(
        target_arch = "x86_64",
        not(target_feature = "avx2"),
        target_feature = "sse2"
    )
)))]
#[inline(always)]
pub(crate) fn xor_phase_with_xor(phase: &mut [u64], left: &[u64], right: &[u64]) {
    debug_assert_eq!(phase.len(), left.len());
    debug_assert_eq!(left.len(), right.len());
    xor_phase_with_xor_scalar(phase, left, right);
}

#[inline(always)]
pub(crate) fn xor_phase_with_scalar(phase: &mut [u64], bits: &[u64]) {
    debug_assert_eq!(phase.len(), bits.len());
    let mut i = 0;
    let len = phase.len();
    while i + 4 <= len {
        phase[i] ^= bits[i];
        phase[i + 1] ^= bits[i + 1];
        phase[i + 2] ^= bits[i + 2];
        phase[i + 3] ^= bits[i + 3];
        i += 4;
    }
    for j in i..len {
        phase[j] ^= bits[j];
    }
}

#[inline(always)]
pub(crate) fn xor_phase_with_xor_scalar(phase: &mut [u64], left: &[u64], right: &[u64]) {
    debug_assert_eq!(phase.len(), left.len());
    debug_assert_eq!(left.len(), right.len());
    let mut i = 0;
    let len = phase.len();
    while i + 4 <= len {
        phase[i] ^= left[i] ^ right[i];
        phase[i + 1] ^= left[i + 1] ^ right[i + 1];
        phase[i + 2] ^= left[i + 2] ^ right[i + 2];
        phase[i + 3] ^= left[i + 3] ^ right[i + 3];
        i += 4;
    }
    for j in i..len {
        phase[j] ^= left[j] ^ right[j];
    }
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
unsafe fn xor_with_avx2(dst: &mut [u64], src: &[u64]) {
    use core::arch::x86_64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 8 <= len {
        let d0 = _mm256_loadu_si256(dst.as_ptr().add(i) as *const __m256i);
        let s0 = _mm256_loadu_si256(src.as_ptr().add(i) as *const __m256i);
        let r0 = _mm256_xor_si256(d0, s0);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i) as *mut __m256i, r0);

        let d1 = _mm256_loadu_si256(dst.as_ptr().add(i + 4) as *const __m256i);
        let s1 = _mm256_loadu_si256(src.as_ptr().add(i + 4) as *const __m256i);
        let r1 = _mm256_xor_si256(d1, s1);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i + 4) as *mut __m256i, r1);
        i += 8;
    }
    while i + 4 <= len {
        let d = _mm256_loadu_si256(dst.as_ptr().add(i) as *const __m256i);
        let s = _mm256_loadu_si256(src.as_ptr().add(i) as *const __m256i);
        let r = _mm256_xor_si256(d, s);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i) as *mut __m256i, r);
        i += 4;
    }
    for j in i..len {
        dst[j] ^= src[j];
    }
}

#[cfg(all(target_arch = "x86_64", target_feature = "avx2"))]
unsafe fn xor_with_xor_avx2(dst: &mut [u64], left: &[u64], right: &[u64]) {
    use core::arch::x86_64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 8 <= len {
        let d0 = _mm256_loadu_si256(dst.as_ptr().add(i) as *const __m256i);
        let l0 = _mm256_loadu_si256(left.as_ptr().add(i) as *const __m256i);
        let r0 = _mm256_loadu_si256(right.as_ptr().add(i) as *const __m256i);
        let t0 = _mm256_xor_si256(l0, r0);
        let out0 = _mm256_xor_si256(d0, t0);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i) as *mut __m256i, out0);

        let d1 = _mm256_loadu_si256(dst.as_ptr().add(i + 4) as *const __m256i);
        let l1 = _mm256_loadu_si256(left.as_ptr().add(i + 4) as *const __m256i);
        let r1 = _mm256_loadu_si256(right.as_ptr().add(i + 4) as *const __m256i);
        let t1 = _mm256_xor_si256(l1, r1);
        let out1 = _mm256_xor_si256(d1, t1);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i + 4) as *mut __m256i, out1);
        i += 8;
    }
    while i + 4 <= len {
        let d = _mm256_loadu_si256(dst.as_ptr().add(i) as *const __m256i);
        let l = _mm256_loadu_si256(left.as_ptr().add(i) as *const __m256i);
        let r = _mm256_loadu_si256(right.as_ptr().add(i) as *const __m256i);
        let t = _mm256_xor_si256(l, r);
        let out = _mm256_xor_si256(d, t);
        _mm256_storeu_si256(dst.as_mut_ptr().add(i) as *mut __m256i, out);
        i += 4;
    }
    for j in i..len {
        dst[j] ^= left[j] ^ right[j];
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[target_feature(enable = "neon")]
unsafe fn xor_with_neon(dst: &mut [u64], src: &[u64]) {
    use core::arch::aarch64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 4 <= len {
        unsafe {
            let d0 = vld1q_u64(dst.as_ptr().add(i));
            let s0 = vld1q_u64(src.as_ptr().add(i));
            let r0 = veorq_u64(d0, s0);
            vst1q_u64(dst.as_mut_ptr().add(i), r0);
            let d1 = vld1q_u64(dst.as_ptr().add(i + 2));
            let s1 = vld1q_u64(src.as_ptr().add(i + 2));
            let r1 = veorq_u64(d1, s1);
            vst1q_u64(dst.as_mut_ptr().add(i + 2), r1);
        }
        i += 4;
    }
    while i + 2 <= len {
        unsafe {
            let d = vld1q_u64(dst.as_ptr().add(i));
            let s = vld1q_u64(src.as_ptr().add(i));
            let r = veorq_u64(d, s);
            vst1q_u64(dst.as_mut_ptr().add(i), r);
        }
        i += 2;
    }
    for j in i..len {
        dst[j] ^= src[j];
    }
}

#[cfg(all(target_arch = "aarch64", target_feature = "neon"))]
#[target_feature(enable = "neon")]
unsafe fn xor_with_xor_neon(dst: &mut [u64], left: &[u64], right: &[u64]) {
    use core::arch::aarch64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 4 <= len {
        unsafe {
            let d0 = vld1q_u64(dst.as_ptr().add(i));
            let l0 = vld1q_u64(left.as_ptr().add(i));
            let r0 = vld1q_u64(right.as_ptr().add(i));
            let t0 = veorq_u64(l0, r0);
            let out0 = veorq_u64(d0, t0);
            vst1q_u64(dst.as_mut_ptr().add(i), out0);
            let d1 = vld1q_u64(dst.as_ptr().add(i + 2));
            let l1 = vld1q_u64(left.as_ptr().add(i + 2));
            let r1 = vld1q_u64(right.as_ptr().add(i + 2));
            let t1 = veorq_u64(l1, r1);
            let out1 = veorq_u64(d1, t1);
            vst1q_u64(dst.as_mut_ptr().add(i + 2), out1);
        }
        i += 4;
    }
    while i + 2 <= len {
        unsafe {
            let d = vld1q_u64(dst.as_ptr().add(i));
            let l = vld1q_u64(left.as_ptr().add(i));
            let r = vld1q_u64(right.as_ptr().add(i));
            let t = veorq_u64(l, r);
            let out = veorq_u64(d, t);
            vst1q_u64(dst.as_mut_ptr().add(i), out);
        }
        i += 2;
    }
    for j in i..len {
        dst[j] ^= left[j] ^ right[j];
    }
}

#[cfg(all(
    target_arch = "x86_64",
    not(target_feature = "avx2"),
    target_feature = "sse2"
))]
unsafe fn xor_with_sse2(dst: &mut [u64], src: &[u64]) {
    use core::arch::x86_64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 4 <= len {
        let d0 = _mm_loadu_si128(dst.as_ptr().add(i) as *const __m128i);
        let s0 = _mm_loadu_si128(src.as_ptr().add(i) as *const __m128i);
        let r0 = _mm_xor_si128(d0, s0);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, r0);

        let d1 = _mm_loadu_si128(dst.as_ptr().add(i + 2) as *const __m128i);
        let s1 = _mm_loadu_si128(src.as_ptr().add(i + 2) as *const __m128i);
        let r1 = _mm_xor_si128(d1, s1);
        _mm_storeu_si128(dst.as_mut_ptr().add(i + 2) as *mut __m128i, r1);
        i += 4;
    }
    while i + 2 <= len {
        let d = _mm_loadu_si128(dst.as_ptr().add(i) as *const __m128i);
        let s = _mm_loadu_si128(src.as_ptr().add(i) as *const __m128i);
        let r = _mm_xor_si128(d, s);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, r);
        i += 2;
    }
    for j in i..len {
        dst[j] ^= src[j];
    }
}

#[cfg(all(
    target_arch = "x86_64",
    not(target_feature = "avx2"),
    target_feature = "sse2"
))]
unsafe fn xor_with_xor_sse2(dst: &mut [u64], left: &[u64], right: &[u64]) {
    use core::arch::x86_64::*;

    let mut i = 0;
    let len = dst.len();
    while i + 4 <= len {
        let d0 = _mm_loadu_si128(dst.as_ptr().add(i) as *const __m128i);
        let l0 = _mm_loadu_si128(left.as_ptr().add(i) as *const __m128i);
        let r0 = _mm_loadu_si128(right.as_ptr().add(i) as *const __m128i);
        let t0 = _mm_xor_si128(l0, r0);
        let out0 = _mm_xor_si128(d0, t0);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, out0);

        let d1 = _mm_loadu_si128(dst.as_ptr().add(i + 2) as *const __m128i);
        let l1 = _mm_loadu_si128(left.as_ptr().add(i + 2) as *const __m128i);
        let r1 = _mm_loadu_si128(right.as_ptr().add(i + 2) as *const __m128i);
        let t1 = _mm_xor_si128(l1, r1);
        let out1 = _mm_xor_si128(d1, t1);
        _mm_storeu_si128(dst.as_mut_ptr().add(i + 2) as *mut __m128i, out1);
        i += 4;
    }
    while i + 2 <= len {
        let d = _mm_loadu_si128(dst.as_ptr().add(i) as *const __m128i);
        let l = _mm_loadu_si128(left.as_ptr().add(i) as *const __m128i);
        let r = _mm_loadu_si128(right.as_ptr().add(i) as *const __m128i);
        let t = _mm_xor_si128(l, r);
        let out = _mm_xor_si128(d, t);
        _mm_storeu_si128(dst.as_mut_ptr().add(i) as *mut __m128i, out);
        i += 2;
    }
    for j in i..len {
        dst[j] ^= left[j] ^ right[j];
    }
}
