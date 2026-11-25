pub mod forward_euler;


fn main() -> Result<(), std::io::Error> {
    forward_euler::eval(0.01, 10.0)
}

// 0.000001