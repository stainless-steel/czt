//! [Chirp Z-transform][1] algorithm.
//!
//! [1]: https://en.wikipedia.org/wiki/Chirp_Z-transform

extern crate complex;
extern crate fft;

use complex::c64;

macro_rules! add_padding(
    ($buffer:expr, $value:expr) => ({
        let buffer = $buffer;
        let capacity = buffer.capacity();
        while buffer.len() < capacity {
            buffer.push($value);
        }
    });
);

macro_rules! increase_to_power_of_two(
    ($number:expr) => ({
        let number = $number;
        assert!(number > 0);
        1 << ((number as f64).log2().ceil() as usize)
    });
);

pub fn forward(data: &[f64]) -> Vec<c64> {
    const ONE: c64 = c64(1.0, 0.0);
    const ZERO: c64 = c64(0.0, 0.0);

    let size = data.len();

    let chirp = {
        use std::f64::consts::PI;
        let theta = -2.0 * PI / size as f64;
        ((-(size as isize) + 1)..(size as isize)).map(|i| {
            let argument = theta * (i * i) as f64 / 2.0;
            c64(argument.cos(), argument.sin())
        }).collect::<Vec<_>>()
    };

    let n = increase_to_power_of_two!(2 * size - 1);

    let mut buffer1 = Vec::with_capacity(n);
    for i in 0..size {
        buffer1.push(chirp[size + i - 1] * data[i]);
    }
    add_padding!(&mut buffer1, ZERO);

    let mut buffer2 = Vec::with_capacity(n);
    for i in 0..(2 * size - 1) {
        buffer2.push(ONE / chirp[i]);
    }
    add_padding!(&mut buffer2, ZERO);

    fft::complex::forward(&mut buffer1);
    fft::complex::forward(&mut buffer2);

    for i in 0..n {
        buffer1[i] = buffer1[i] * buffer2[i];
    }

    fft::complex::inverse(&mut buffer1);

    ((size - 1)..(2 * size - 1)).map(|i| buffer1[i] * chirp[i]).collect()
}
