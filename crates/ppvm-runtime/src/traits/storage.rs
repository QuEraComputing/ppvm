use bitvec::view::BitViewSized;
use std::hash::Hash;

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
