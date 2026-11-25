use std::fs::File;
use std::io::{self, BufWriter, Write};

 /// Implementation of the **Forward Euler** method *(explizites Euler-Verfahren)* as a **first order** solver.
 /// 
 /// This uses default values for now as well as a simplified unit circle approach. No masses or constants.
 /// 
 /// ### Defaults
 /// * r_0 = (1.0, 0.0, 0.0)
 /// * v_0 = (0.0, 0.5, 0.0)
pub fn eval(timestep: f64, duration: f64) -> io::Result<()> {
    let dt = timestep;
    let t_max: u64 = f64::ceil(duration / dt) as u64;

    let mut r: [f64; 3] = [1.0, 0.0, 0.0];
    let mut v: [f64; 3] = [0.0, 0.5, 0.0];

    let file = File::create("output.txt")?;
    let mut stream = BufWriter::new(file);
    writeln!(stream, "#t r_0 r_1 r_2 v_0 v_1 v_2")?;
    writeln!(stream, "0.0 {} {} {} {} {} {}", r[0], r[1], r[2], v[0], v[1], v[2])?;

    let dt_thresh = 1000;
    let mut t_thresh = dt_thresh;
    
    for ts in 0..t_max {
        let t = ts as f64 * dt;
        if ts >= t_thresh {
            println!("{} / {}", t, t_max);
            writeln!(stream, "{} {} {} {} {} {} {}", t, r[0], r[1], r[2], v[0], v[1], v[2])?;
            t_thresh += dt_thresh;
        }
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

        let mut a: [f64; 3] = [0.0; 3];

        for k in 0..3 {
            a[k] = -r[k] / (r2 * r2.sqrt()); 
        }
        for k in 0..3 {
            r[k] += v[k] * dt;
            v[k] += a[k] * dt;
        }
        
    }
    stream.flush()?;
    Ok(())
}