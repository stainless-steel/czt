use num_complex::Complex;
use num_traits::{Float, FloatConst};
use std::ops::Mul;

/// The transform.
pub trait Transform<T> {
    /// Perform the transform.
    ///
    /// ## References
    ///
    /// 1. https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm
    fn transform(&self, m: usize, w: Complex<T>, a: Complex<T>) -> Vec<Complex<T>>;
}

impl<D, T> Transform<T> for [D]
    where D: Copy + Mul<Complex<T>, Output=Complex<T>>, T: Float + FloatConst
{
    fn transform(&self, m: usize, w: Complex<T>, a: Complex<T>) -> Vec<Complex<T>> {
        use dft::{Operation, Plan, Transform};
        use num_traits::{One, Zero};

        let zero = Complex::zero();
        let one = Complex::one();
        let two = T::one() + T::one();
        let n = self.len();
        let (modulus, argument) = w.to_polar();
        let factor = ((-(n as isize) + 1)..(if n > m { n } else { m } as isize)).map(|i| {
            let power = T::from(i * i).unwrap() / two;
            Complex::from_polar(&modulus.powf(power), &(argument * power))
        }).collect::<Vec<_>>();
        let p = (n + m - 1).next_power_of_two();
        let mut buffer1 = vec![zero; p];
        let (modulus, argument) = a.to_polar();
        for i in 0..n {
            let power = -T::from(i).unwrap();
            let a = Complex::from_polar(&modulus.powf(power), &(argument * power));
            buffer1[i] = self[i] * factor[n + i - 1] * a;
        }
        let mut buffer2 = vec![zero; p];
        for i in 0..(n + m - 1) {
            buffer2[i] = one / factor[i];
        }
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
