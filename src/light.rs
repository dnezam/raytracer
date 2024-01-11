use crate::{errors::LightError, Color, Tuple};

type Result<T> = std::result::Result<T, LightError>;

/// Represents a light source with no size (i.e. a point).
#[derive(Debug, Copy, Clone)]
pub struct PointLight {
    // Invariant: Must be a point
    position: Tuple,
    pub intensity: Color,
}

impl PointLight {
    /// Create a point light.
    ///
    /// # Returns
    /// Returns Ok(PointLight) if `position` is a point.
    /// Returns Err(NotAPoint) if `position` is not a point.
    pub fn new(position: Tuple, intensity: Color) -> Result<Self> {
        if !position.is_point() {
            Err(LightError::NotAPoint)
        } else {
            Ok(Self {
                position,
                intensity,
            })
        }
    }

    /// Returns the position.
    pub fn position(self) -> Tuple {
        self.position
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let intensity = Color::new(1.0, 1.0, 1.0);
        let position = Tuple::point(0.0, 0.0, 0.0);
        let light = PointLight::new(position, intensity).unwrap();

        assert_eq!(light.position(), position);
        assert_eq!(light.intensity, intensity);
    }
}
