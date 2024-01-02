use raytracer::Tuple;

fn main() {
    let mut p = Projectile {
        position: Tuple::point(0.0, 1.0, 0.0),
        velocity: Tuple::vector(1.0, 1.0, 0.0).normalize(),
    };
    let e = Environment {
        gravity: Tuple::vector(0.0, -0.9, 0.0),
        wind: Tuple::vector(-0.01, 0.0, 0.0),
    };

    while p.position.y > 0.0 {
        println!("{:?}", p);
        p = tick(e, p);
    }
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
