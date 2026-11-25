#[cfg(feature = "single-precision")]
pub type Real = f32;

#[cfg(not(feature = "single-precision"))]
pub type Real = f64;

pub const G: Real = 6.674_30e-11;

/// Implementing a Vector and related methods needed for the solvers.
/// 
/// Currently only for forward_euler.rs. No mass.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vector{
    pub x: Real, 
    pub y: Real, 
    pub z: Real,
}

impl Vector{
    pub fn new(x: Real, y: Real, z: Real) -> Self {
        Self { x, y, z }
    }

    pub fn mag_sq(&self) -> Real {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn mag(&self) -> Real {
        self.mag_sq().sqrt()
    }

    pub fn kinetic_energy(&self) -> Real{
        0.5 * self.mag_sq()
    } 

    pub fn potential_energy(&self) -> Real{
        -1.0 / self.mag_sq().sqrt()
    }

    pub fn total_energy(v: &Vector, r: &Vector) -> Real {
        v.kinetic_energy() + r.potential_energy()
    }

    pub fn new_acceleration(distance: &Vector) -> Self {
        let r_2 = distance.mag_sq();
        let mag = r_2.sqrt();
        Self {
            x: -distance.x / (r_2 * mag),
            y: -distance.y / (r_2 * mag),
            z: -distance.z / (r_2 * mag),
        }
    }

    /// Calculates a timestep over the vector for Forward Euler.
    /// ### Arguments
    /// * integral - takes (a) **v** if &self == **r** or (b) **a** if &self == **v**
    /// * dt - timestep length in seconds
    pub fn forward_euler_step(&mut self, derivative: &Vector, dt: &Real) {
        self.x += derivative.x * dt;
        self.y += derivative.y * dt;
        self.z += derivative.z * dt;
    }
}