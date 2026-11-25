use std::fs::File;
use std::io::{self, BufWriter, Write};
use orbit_simulation::{Real, Vector};

/// Implementation of the **Forward Euler** method *(explizites Euler-Verfahren)* as a **first order** solver.
/// 
/// This uses default values for now as well as a simplified unit circle approach. No masses or constants.
/// 
/// ### Defaults
/// * r_0 = (1.0, 0.0, 0.0)
/// * v_0 = (0.0, 0.5, 0.0)
pub fn eval(dt: Real, duration: Real) -> io::Result<()> {
    let t_max: u64 = Real::ceil(duration / dt) as u64;

    let mut r = Vector::new(
        1.0, 
        0.0, 
        0.0,
    );
    let mut v = Vector::new(
        0.0, 
        0.5, 
        0.0,
    );

    let file = File::create("output.csv")?;
    let mut stream = BufWriter::new(file);

    let e_in = Vector::total_energy(&v, &r);
    println!("Initial total energy is {}", e_in);

    writeln!(stream, "#t r_0 r_1 r_2 v_0 v_1 v_2 total_energy")?;
    writeln!(stream, "0.0 {} {} {} {} {} {} {}", r.x, r.y, r.z, v.x, v.y, v.z, e_in)?;


    let dt_thresh = t_max / 1000;
    let mut t_thresh = dt_thresh;
    
    for ts in 0..t_max {
        let t = ts as Real * dt;
        if ts >= t_thresh {
            let e = Vector::total_energy(&v, &r);
            // println!("{} / {}", t, t_max);
            writeln!(stream, "{} {} {} {} {} {} {} {}", t, r.x, r.y, r.z, v.x, v.y, v.z, e)?;
            t_thresh += dt_thresh;
        }

        let a= Vector::new_acceleration(&r);
        r.forward_euler_step(&v, &dt);
        v.forward_euler_step(&a, &dt);
    }

    let e = Vector::total_energy(&v, &r);

    println!("Final total energy is {}", e);
    println!("absolute energy error : {}", e-e_in);
    println!("relative energy error : {}", (e - e_in) / e_in);

    stream.flush()?;
    Ok(())
}