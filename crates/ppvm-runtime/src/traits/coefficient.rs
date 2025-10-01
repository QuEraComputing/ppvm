use num::traits::{Float, One, Zero};
use std::fmt::{Debug, Display};

pub trait Coefficient:
    Clone
    + Display
    + Debug
    + Send
    + Sync
    + One
    + Zero
    + PartialEq
    + std::ops::Neg<Output = Self>
    // + std::ops::Add
    // + std::ops::Sub
    // + std::ops::Mul
    + std::ops::AddAssign
    + std::ops::SubAssign
    + std::ops::MulAssign
{
    fn mul_sign(&self, sign: i8) -> Self;
    fn half(&self) -> Self;
    fn sin_cos(&self) -> (Self, Self);
}

impl Coefficient for f64 {
    fn mul_sign(&self, sign: i8) -> Self {
        (sign as f64) * (*self)
    }

    fn half(&self) -> Self {
        *self / 2.0
    }

    fn sin_cos(&self) -> (Self, Self) {
        Float::sin_cos(*self)
    }
}