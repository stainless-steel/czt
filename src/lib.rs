//! [Chirp Z-transform][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Chirp_Z-transform

extern crate dft;
extern crate num_complex;
extern crate num_traits;

use num_complex::Complex;

mod transform;

pub use transform::Transform;

/// A complex number with 32-bit parts.
#[allow(non_camel_case_types)]
pub type c32 = Complex<f32>;

/// A complex number with 64-bit parts.
#[allow(non_camel_case_types)]
pub type c64 = Complex<f64>;

/// Perform the transform.
///
/// The function is a shortcut for `Transform::transform`.
#[inline(always)]
pub fn transform<D: ?Sized, T>(data: &D, m: usize, w: Complex<T>, a: Complex<T>) -> Vec<Complex<T>>
    where D: Transform<T>
{
    Transform::transform(data, m, w, a)
}
