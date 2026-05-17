use bitvec::view::BitViewSized;
use std::hash::Hash;

/// Backing storage for a [`PauliWord`](crate::word::PauliWord) — a
/// fixed-size, `Copy`-able block of bits (typically `[u8; N]` or
/// `[u64; N]`).
pub trait PauliStorage:
    BitViewSized
    + Clone
    + Copy
    + Hash
    + Eq
    + PartialEq
    + Send
    + Sync
    + std::fmt::Debug
    + Ord
    + PartialOrd
{
}
impl<
    A: BitViewSized
        + Clone
        + Copy
        + Hash
        + Eq
        + PartialEq
        + Send
        + Sync
        + std::fmt::Debug
        + Ord
        + PartialOrd,
> PauliStorage for A
{
}
