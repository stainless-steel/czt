//! [Chirp Z-transform][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Chirp_Z-transform

extern crate dft;
extern crate num_complex;
extern crate num_traits;

use num_complex::Complex;
use num_traits::Float;
use std::ops::Mul;

/// A complex number with 32-bit parts.
#[allow(non_camel_case_types)]
pub type c32 = Complex<f32>;

/// A complex number with 64-bit parts.
#[allow(non_camel_case_types)]
pub type c64 = Complex<f64>;

/// The transform.
pub trait Transform<T> {
    /// Perform the transform.
    ///
    /// ## References
    ///
    /// 1. https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm
    fn transform(&self, m: usize, w: Complex<T>, a: Complex<T>) -> Vec<Complex<T>>;
}

macro_rules! add_padding(
    ($buffer:expr) => ({
        let buffer = $buffer;
        let zero = Complex::zero();
        let capacity = buffer.capacity();
        while buffer.len() < capacity {
            buffer.push(zero);
        }
    });
);

impl<D, T> Transform<T> for [D] where D: Copy + Mul<Complex<T>, Output=Complex<T>>, T: Float {
    fn transform(&self, m: usize, w: Complex<T>, a: Complex<T>) -> Vec<Complex<T>> {
        use dft::{Operation, Plan, Transform};
        use num_traits::{One, Zero};

        let n = self.len();
        let factor = {
            let two = T::one() + T::one();
            let (modulus, argument) = w.to_polar();
            ((-(n as isize) + 1)..(if n > m { n } else { m } as isize)).map(|i| {
                let power = T::from(i * i).unwrap() / two;
                Complex::from_polar(&modulus.powf(power), &(argument * power))
            }).collect::<Vec<_>>()
        };
        let p = (n + m - 1).next_power_of_two();
        let mut buffer1 = Vec::with_capacity(p);
        {
            let (modulus, argument) = a.to_polar();
            for i in 0..n {
                let power = -T::from(i).unwrap();
                let a = Complex::from_polar(&modulus.powf(power), &(argument * power));
                buffer1.push(self[i] * factor[n + i - 1] * a);
            }
        }
        add_padding!(&mut buffer1);
        let one = Complex::one();
        let mut buffer2 = Vec::with_capacity(p);
        for i in 0..(n + m - 1) {
            buffer2.push(one / factor[i]);
        }
        add_padding!(&mut buffer2);
        let plan = Plan::new(Operation::Forward, p);
        buffer1.transform(&plan);
        buffer2.transform(&plan);
        for i in 0..p {
            buffer1[i] = buffer1[i] * buffer2[i];
        }
        let plan = Plan::new(Operation::Inverse, p);
        buffer1.transform(&plan);
        ((n - 1)..(n + m - 1)).map(|i| buffer1[i] * factor[i]).collect()
    }
}
