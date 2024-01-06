use std::path::Path;

use raytracer::{Canvas, Color, Tuple};

fn main() {
    let width = 950;
    let height = 550;
    let mut c = Canvas::new(width, height);
    let mut p = Projectile {
        position: Tuple::point(0.0, 0.5, 0.0),
        // Normalizing a vector returns a vector.
        velocity: Tuple::vector(1.0, 1.8, 0.0).normalize().unwrap() * 11.5,
    };
    let e = Environment {
        gravity: Tuple::vector(0.0, -0.1, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    while p.position.y > 0.0 {
        let x = p.position.x.round() as usize;
        let y = p.position.y.round() as usize;

        let Ok(canvas_position) = c.get_mut(x, height - y) else {
            continue;
        };
        *canvas_position = Color::new(1.0, 0.0, 0.0);

        p = tick(e, p);
    }

    c.save(Path::new("./projectile.ppm"))
        .expect("Unable to save canvas.")
}

#[derive(Debug, Copy, Clone)]
struct Projectile {
    position: Tuple, // point
    velocity: Tuple, // vector
}

#[derive(Debug, Copy, Clone)]
struct Environment {
    gravity: Tuple, // vector
    wind: Tuple,    // vector
}

fn tick(env: Environment, proj: Projectile) -> Projectile {
    let position = proj.position + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    Projectile { position, velocity }
}
