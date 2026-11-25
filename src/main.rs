pub mod forward_euler;

pub const G: f64 = 6.674_30e-11;

fn main() {
    let _ = forward_euler::eval(0.000001, 10.0);
}
