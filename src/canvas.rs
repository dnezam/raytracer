use std::ops::{Index, IndexMut};

use crate::color::Color;

#[derive(PartialEq)]
/// A rectangular grid of pixels.
pub struct Canvas {
    width: usize,
    height: usize,
    pixels: Vec<Color>,
}

impl Index<[usize; 2]> for Canvas {
    type Output = Color;

    /// Implements indexing of the form canvas[[x, y]] in immutable contexts.
    ///
    /// The x and y parameters are assumed to be 0-based: x may be anywhere from
    /// 0 to width - 1 (inclusive), and y may be anywhere from 0 to height - 1 (inclusive).
    ///
    /// # Panics
    /// `index` will panic if x or y is out-of-bounds.
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(index[0] < self.width);
        assert!(index[1] < self.height);
        &self.pixels[index[1] * self.width + index[0]]
    }
}

impl IndexMut<[usize; 2]> for Canvas {
    /// Implements indexing of the form canvas[[x, y]] in mutable contexts.
    ///
    /// The x and y parameters are assumed to be 0-based: x may be anywhere from
    /// 0 to width - 1 (inclusive), and y may be anywhere from 0 to height - 1 (inclusive).
    ///
    /// # Panics
    /// `index` will panic if x or y is out-of-bounds.
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(index[0] < self.width);
        assert!(index[1] < self.height);
        &mut self.pixels[index[1] * self.width + index[0]]
    }
}

impl Canvas {
    /// Create a new Canvas.
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            pixels: vec![Color::new(0.0, 0.0, 0.0); width * height],
        }
    }

    /// Returns an iterator over the Canvas.
    pub fn iter(&self) -> impl Iterator<Item = &Color> {
        self.pixels.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width, 10);
        assert_eq!(c.height, 20);
        for point in c.iter() {
            assert_eq!(*point, Color::new(0.0, 0.0, 0.0))
        }
    }

    #[test]
    fn write() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);
        c[[2, 3]] = red;
        assert_eq!(c[[2, 3]], red);
    }
}
