use std::ops::{Add, Mul, Sub};

use crate::utils;

/// Used to represent color in the RGB format.
///
/// This is separate from Color, as Color stores values between 0 and 1, while Rgb
/// stores values between 0 to 255.
#[derive(Debug, Copy, Clone)]
pub struct Rgb(pub u8, pub u8, pub u8);

/// Used to represent color.
///
/// Usually, values are between 0 and 1, with 0 meaning the color is entirely
/// absent, and 1 meaning the color is fully present. However, we don't constrain
/// them immediately, as the value may leave this interval while going through
/// transformations. If we were to limit the color prematurely, we run this risk
/// of making parts of our scene to bright or dark in the final image.
#[derive(Debug, Copy, Clone)]
pub struct Color {
    /// The red component.
    pub red: f64,
    /// The green component.
    pub green: f64,
    /// The blue component.
    pub blue: f64,
}

impl PartialEq for Color {
    /// Implements approximate equality for colors.
    fn eq(&self, other: &Self) -> bool {
        utils::eq(self.red, other.red)
            && utils::eq(self.green, other.green)
            && utils::eq(self.blue, other.blue)
    }
}

impl Add for Color {
    type Output = Self;

    /// Implements the addition of two colors.
    fn add(self, other: Self) -> Self::Output {
        Self {
            red: self.red + other.red,
            green: self.green + other.green,
            blue: self.blue + other.blue,
        }
    }
}

impl Sub for Color {
    type Output = Self;

    /// Implements the subtraction of two colors.
    fn sub(self, other: Self) -> Self::Output {
        Self {
            red: self.red - other.red,
            green: self.green - other.green,
            blue: self.blue - other.blue,
        }
    }
}

/// Scalar multiplication of the form Color * f64.
impl Mul<f64> for Color {
    type Output = Self;

    /// Implements scalar multiplication of the form Color * f64.
    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            red: self.red * rhs,
            green: self.green * rhs,
            blue: self.blue * rhs,
        }
    }
}

/// Scalar multiplication of the form f64 * Color.
impl Mul<Color> for f64 {
    type Output = Color;

    /// Implements scalar multiplication of the form f64 * Color.
    fn mul(self, rhs: Color) -> Self::Output {
        Self::Output {
            red: self * rhs.red,
            green: self * rhs.green,
            blue: self * rhs.blue,
        }
    }
}

// Multiplication of the form Color * Color.
impl Mul<Color> for Color {
    type Output = Self;

    /// Implements multiplication of the form Color * Color.
    ///
    /// This is technically called the Hadamard product, Schur product
    /// or element-wise product. It is used to blend two colors together.
    fn mul(self, rhs: Color) -> Self::Output {
        Self::Output {
            red: self.red * rhs.red,
            green: self.green * rhs.green,
            blue: self.blue * rhs.blue,
        }
    }
}

impl Color {
    /// Create a new color.
    pub fn new(red: f64, green: f64, blue: f64) -> Color {
        Color { red, green, blue }
    }

    /// Converts Color into Rgb.
    pub fn to_rgb(self) -> Rgb {
        Rgb(
            (self.red * 255.0).round() as u8,
            (self.green * 255.0).round() as u8,
            (self.blue * 255.0).round() as u8,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let c = Color::new(-0.5, 0.4, 1.7);
        assert!(utils::eq(c.red, -0.5));
        assert!(utils::eq(c.green, 0.4));
        assert!(utils::eq(c.blue, 1.7));
    }

    #[test]
    fn add() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 + c2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn subtract() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(c1 - c2, Color::new(0.2, 0.5, 0.5));
    }

    #[test]
    fn multiply_scalar() {
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c * 2.0, Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn multiply_scalar_commutative() {
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c * 2.0, 2.0 * c);
    }

    #[test]
    fn multiply_colors() {
        let c1 = Color::new(1.0, 0.2, 0.4);
        let c2 = Color::new(0.9, 1.0, 0.1);
        assert_eq!(c1 * c2, Color::new(0.9, 0.2, 0.04));
    }
}
