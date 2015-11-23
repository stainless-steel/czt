//! [Chirp Z-transform][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Chirp_Z-transform

extern crate dft;

pub use dft::{Operation, Plan, Transform, c64};
use std::ops::Mul;

macro_rules! add_padding(
    ($buffer:expr, $value:expr) => ({
        let buffer = $buffer;
        let capacity = buffer.capacity();
        while buffer.len() < capacity {
            buffer.push($value);
        }
    });
);

/// Perform the forward transformation.
///
/// ## References
///
/// 1. https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm
pub fn forward<T>(data: &[T], m: usize, w: c64, a: c64) -> Vec<c64>
    where T: Mul<c64, Output=c64> + Copy
{
    const ONE: c64 = c64 { re: 1.0, im: 0.0 };
    const ZERO: c64 = c64 { re: 0.0, im: 0.0 };

    let n = data.len();

    let factor = {
        let (modulus, argument) = w.to_polar();
        ((-(n as isize) + 1)..(if n > m { n } else { m } as isize)).map(|i| {
            let power = (i * i) as f64 / 2.0;
            c64::from_polar(&modulus.powf(power), &(argument * power))
        }).collect::<Vec<_>>()
    };

    let p = (n + m - 1).next_power_of_two();

    let mut buffer1 = Vec::with_capacity(p);
    {
        let (modulus, argument) = a.to_polar();
        for i in 0..n {
            let power = -(i as f64);
            let a = c64::from_polar(&modulus.powf(power), &(argument * power));
            buffer1.push(data[i] * factor[n + i - 1] * a);
        }
    }
    add_padding!(&mut buffer1, ZERO);

    let mut buffer2 = Vec::with_capacity(p);
    for i in 0..(n + m - 1) {
        buffer2.push(ONE / factor[i]);
    }
    add_padding!(&mut buffer2, ZERO);

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
