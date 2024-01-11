use crate::{errors::MaterialError, utils, Color};

type Result<T> = std::result::Result<T, MaterialError>;

/// Encapsulates color and the attributes from the Phong reflection model.
#[derive(Debug, Copy, Clone)]
pub struct Material {
    pub color: Color,
    // Invariant: Must be nonnegative
    ambient: f64,
    // Invariant: Must be nonnegative
    diffuse: f64,
    // Invariant: Must be nonnegative
    specular: f64,
    // Invariant: Must be nonnegative
    shininess: f64,
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.color == other.color
            && utils::eq(self.ambient, other.ambient)
            && utils::eq(self.diffuse, other.diffuse)
            && utils::eq(self.specular, other.specular)
            && utils::eq(self.shininess, other.shininess)
    }
}

impl Default for Material {
    fn default() -> Self {
        Self {
            color: Color::new(1.0, 1.0, 1.0),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
        }
    }
}

impl Material {
    /// Sets the value for `ambient`.
    ///
    /// Typical values are between 0 and 1.
    ///
    /// # Returns
    /// Returns Ok(()) if `ambient` is a nonnegative float.
    /// Returns Err(MaterialError::NegativeFloat) if `ambient` is a negative float.
    pub fn set_ambient(&mut self, ambient: f64) -> Result<()> {
        if ambient < 0.0 {
            Err(MaterialError::NegativeFloat)
        } else {
            Ok(self.ambient = ambient)
        }
    }

    /// Sets the value for `diffuse`.
    ///
    /// Typical values are between 0 and 1.
    ///
    /// # Returns
    /// Returns Ok(()) if `diffuse` is a nonnegative float.
    /// Returns Err(MaterialError::NegativeFloat) if `diffuse` is a negative float.
    pub fn set_diffuse(&mut self, diffuse: f64) -> Result<()> {
        if diffuse < 0.0 {
            Err(MaterialError::NegativeFloat)
        } else {
            Ok(self.diffuse = diffuse)
        }
    }

    /// Sets the value for `specular`.
    ///
    /// Typical values are between 0 and 1.
    ///
    /// # Returns
    /// Returns Ok(()) if `specular` is a nonnegative float.
    /// Returns Err(MaterialError::NegativeFloat) if `specular` is a negative float.
    pub fn set_specular(&mut self, specular: f64) -> Result<()> {
        if specular < 0.0 {
            Err(MaterialError::NegativeFloat)
        } else {
            Ok(self.specular = specular)
        }
    }

    /// Sets the value for `shininess`.
    ///
    /// Typical values are between 10 (very large highlight) and 200 (very small highlight).
    ///
    /// # Returns
    /// Returns Ok(()) if `shininess` is a nonnegative float.
    /// Returns Err(MaterialError::NegativeFloat) if `shininess` is a negative float.
    pub fn set_shininess(&mut self, shininess: f64) -> Result<()> {
        if shininess < 0.0 {
            Err(MaterialError::NegativeFloat)
        } else {
            Ok(self.shininess = shininess)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default() {
        let m = Material::default();
        assert_eq!(m.color, Color::new(1.0, 1.0, 1.0));
        assert!(utils::eq(m.ambient, 0.1));
        assert!(utils::eq(m.diffuse, 0.9));
        assert!(utils::eq(m.specular, 0.9));
        assert!(utils::eq(m.shininess, 200.0));
    }
}
