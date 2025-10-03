use num::traits::Float;

pub trait Coefficient:
    PartialEq
    + Clone
    + num::Zero
    + From<i32>
    + From<f32>
    + From<f64>
    + std::ops::Neg<Output = Self>
    + std::ops::Add<f64, Output = Self>
    + std::ops::Add<Self, Output = Self>
    + std::ops::Sub<f64, Output = Self>
    + std::ops::Sub<Self, Output = Self>
    + std::ops::Mul<f64, Output = Self>
    + std::ops::Mul<Self, Output = Self>
    + std::ops::AddAssign<f64>
    + std::ops::AddAssign<Self>
    + std::ops::MulAssign<f64>
    + std::ops::MulAssign<Self>
    + std::iter::Sum
    + Sync
    + Send
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
