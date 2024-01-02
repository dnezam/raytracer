// TODO Check if we can loops more idiomatic/"abstract"
use std::{
    ops::{Index, IndexMut},
    slice::Chunks,
};

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
    /// Create a new canvas.
    pub fn new(width: usize, height: usize) -> Canvas {
        Canvas {
            width,
            height,
            pixels: vec![Color::new(0.0, 0.0, 0.0); width * height],
        }
    }

    /// Returns an iterator over all pixels.
    pub fn iter(&self) -> impl Iterator<Item = &Color> {
        self.pixels.iter()
    }

    /// Returns an iterator over the rows.
    pub fn rows(&self) -> Chunks<Color> {
        self.pixels.chunks(self.width)
    }

    /// Returns a mutable iterator over all pixels.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Color> {
        self.pixels.iter_mut()
    }

    /// Converts the canvas into the PPM format.
    fn to_ppm(&self) -> String {
        let header = format!("P3\n{} {}\n255", self.width, self.height);
        let data = self
            .rows()
            // .., Color {1.0, 0.0, 0.0}, .. -> .., ["255", "0", "0"]
            .map(|row| {
                row.iter().flat_map(|color| {
                    let rgb = color.to_rgb();
                    [rgb.0.to_string(), rgb.1.to_string(), rgb.2.to_string()]
                })
            })
            // Combine the numbers in a row to a string with lines of length at most 70
            .map(|row| {
                const MAX_LENGTH: usize = 70;
                let mut result = Vec::new();
                let mut current_line = String::new();

                for number in row {
                    if current_line.is_empty() {
                        assert!(number.len() < MAX_LENGTH);
                        current_line.push_str(&number);
                    } else if current_line.len() + number.len() + 2 <= 70 {
                        // We can only add a number if the current length + space +
                        // the number + \n are at most MAX_LENGTH characters.
                        current_line.push(' ');
                        current_line.push_str(&number);
                    } else {
                        result.push(current_line.clone());
                        current_line.clear();
                        assert!(number.len() < MAX_LENGTH);
                        current_line.push_str(&number);
                    }
                }

                if !current_line.is_empty() {
                    result.push(current_line);
                }

                result.join("\n")
            })
            .collect::<Vec<String>>()
            .join("\n");

        header + "\n" + &data + "\n"
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
        for pixel in c.iter() {
            assert_eq!(*pixel, Color::new(0.0, 0.0, 0.0));
        }
    }

    #[test]
    fn write() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1.0, 0.0, 0.0);
        c[[2, 3]] = red;
        assert_eq!(c[[2, 3]], red);
    }

    #[test]
    fn to_ppm_header() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();

        let expected_lines = ["P3", "5 3", "255"];
        let actual_lines: Vec<&str> = ppm.lines().take(3).collect();
        assert_eq!(actual_lines, expected_lines);
    }

    #[test]
    fn to_ppm_data() {
        let mut c = Canvas::new(5, 3);
        c[[0, 0]] = Color::new(1.5, 0.0, 0.0);
        c[[2, 1]] = Color::new(0.0, 0.5, 0.0);
        c[[4, 2]] = Color::new(-0.5, 0.0, 1.0);
        let ppm = c.to_ppm();

        let expected_lines = [
            "255 0 0 0 0 0 0 0 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 128 0 0 0 0 0 0 0",
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 255",
        ];
        let actual_lines: Vec<&str> = ppm.lines().skip(3).take(3).collect();
        assert_eq!(actual_lines, expected_lines);
    }

    #[test]
    fn to_ppm_split_long_lines() {
        let mut c = Canvas::new(10, 2);
        for pixel in c.iter_mut() {
            *pixel = Color::new(1.0, 0.8, 0.6);
        }
        let ppm = c.to_ppm();

        let expected_lines = [
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204",
            "153 255 204 153 255 204 153 255 204 153 255 204 153",
            "255 204 153 255 204 153 255 204 153 255 204 153 255 204 153 255 204",
            "153 255 204 153 255 204 153 255 204 153 255 204 153",
        ];
        let actual_lines: Vec<&str> = ppm.lines().skip(3).take(4).collect();
        assert_eq!(actual_lines, expected_lines);
    }

    #[test]
    fn to_ppm_terminates_with_newline() {
        let c = Canvas::new(5, 3);
        let ppm = c.to_ppm();
        assert!(ppm.ends_with('\n'))
    }
}
