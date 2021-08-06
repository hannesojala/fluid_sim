// Hate this code so much but IDK how to make it prettier.
macro_rules! i(
    ($n:expr, $x:expr, $y:expr) => {
        ($n * $y + $x) as usize
    }
);

fn fast_get(x: i32, y: i32, n: i32, a: &mut [f32]) -> f32 {
    let i = i!(n,x,y);
    if let Some(v) = a.get(i) {
        *v
    } else {
        0.0
    }
}

fn lerp(v1: f32, v2: f32, k: f32) -> f32 { v1 + k * (v2 - v1) }

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

pub struct Fluid {
    visc: f32,
    diff: f32,
    vort: f32,
    qual: u32,
    size: i32,
    vx: Vec<f32>,
    vy: Vec<f32>,
    dye_r: Vec<f32>,
    dye_g: Vec<f32>,
    dye_b: Vec<f32>,
}

impl Fluid {
    pub fn new(visc: f32, diff: f32, vort: f32, qual: u32, size: i32) -> Fluid {
        let flat_size = ((size+2)*(size+2)) as usize;
        Fluid {
            visc, diff, size, vort, qual,
            vx: vec![0.0; flat_size],
            vy: vec![0.0; flat_size],
            dye_r: vec![0.0; flat_size],
            dye_g: vec![0.0; flat_size],
            dye_b: vec![0.0; flat_size],
        }
    }

    pub fn update(&mut self, dt_s: f32) {
        self.vx = self.advect_field(&self.vx, dt_s);
        self.vy = self.advect_field(&self.vy, dt_s);

        if self.diff > 0. {
            let velocity_diffusion_rate = dt_s * self.visc * ((self.size - 2)*(self.size - 2)) as f32;
            self.vx = self.diffuse_field(&self.vx, velocity_diffusion_rate, true);
            self.vy = self.diffuse_field(&self.vy, velocity_diffusion_rate, true); // bounds
        }
        
        self.confine_vorticity(dt_s);
        self.remove_divergence();

        self.dye_r = self.advect_field(&self.dye_r, dt_s);
        self.dye_g = self.advect_field(&self.dye_g, dt_s);
        self.dye_b = self.advect_field(&self.dye_b, dt_s);
        
        if self.diff > 0. {
            let dye_diffusion_rate = dt_s * self.diff * ((self.size - 2)*(self.size - 2)) as f32;
            self.dye_r = self.diffuse_field(&self.dye_r, dye_diffusion_rate, false);
            self.dye_g = self.diffuse_field(&self.dye_g, dye_diffusion_rate, false);
            self.dye_b = self.diffuse_field(&self.dye_b, dye_diffusion_rate, false);  // bounds
        }
    }

    pub fn set_dye(&mut self, x: i32, y: i32, rgb: (f32, f32, f32)) {
        let i = i!(self.size, x,y);
        self.dye_r[i] = rgb.0;
        self.dye_g[i] = rgb.1;
        self.dye_b[i] = rgb.2;
    }

    pub fn add_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        let i = i!(self.size, x,y);
        self.vx[i] += vx;
        self.vy[i] += vy;
    }

    pub fn dye_at(&self, x: i32, y: i32) -> (f32, f32, f32) {
        let i = i!(self.size, x,y);
        (self.dye_r[i], self.dye_g[i], self.dye_b[i])
    }

    pub fn vel_at(&self, x: i32, y: i32) -> (f32,f32) {
        let i = i!(self.size, x,y);
        (self.vx[i], self.vy[i])
    }

    fn remove_divergence(&mut self) {
        let n = self.size;
        let mut div = vec![0.0; self.vx.len()];
        let mut p = vec![0.0; self.vx.len()];
        for y in 1..n-1 {
            for x in 1..n-1 {
                div[i!(n, x,y)] = -0.5 * (self.vx[i!(n, x+1,y)] - self.vx[i!(n, x-1,y)] + self.vy[i!(n, x,y+1)] - self.vy[i!(n, x,y-1)]) / n as f32;
            }
        }
        Fluid::bound_field(&mut div, n, false);
        for _ in 0..self.qual {
            let p_0 = p.clone(); // This clone actually makes it faster?
            for y in 1..n-1 {
                for x in 1..n-1 {
                    p[i!(n, x,y)] = (div[i!(n, x,y)] + p_0[i!(n, x-1,y)] + p_0[i!(n, x+1,y)] + p_0[i!(n, x,y-1)] + p_0[i!(n, x,y+1)]) / 4.0;
                }
            }
            Fluid::bound_field(&mut p, n, false);
        }
        for y in 1..n-1 {
            for x in 1..n-1 {
                self.vx[i!(n, x,y)] -= 0.5 * (p[i!(n, x+1,y)] - p[i!(n, x-1,y)]) * n as f32;
                self.vy[i!(n, x,y)] -= 0.5 * (p[i!(n, x,y+1)] - p[i!(n, x,y-1)]) * n as f32;
            }
        }
        Fluid::bound_field(&mut self.vx, n, true);
        Fluid::bound_field(&mut self.vy, n, true);
    }

    // Reintroduce vorticity lost by numerical inaccuracy
    fn confine_vorticity(&mut self, dt_s: f32) {
        for y in 2..self.size-2 {
            for x in 2..self.size-2 {
                // Get the gradient of the magnitude of velocity curl
                let mut vort_grad_x = self.vel_curl_at(x, y-1).abs() - self.vel_curl_at(x, y+1).abs();
                let mut vort_grad_y = self.vel_curl_at(x+1, y).abs() - self.vel_curl_at(x-1, y).abs();
                // Normalized and scaled by vorticity constant
                let len = (vort_grad_x*vort_grad_x + vort_grad_y*vort_grad_y).sqrt().max(f32::EPSILON); // max prevents divide by zero
                vort_grad_x *= self.vort / len;
                vort_grad_y *= self.vort / len;
                // Adjust the velocity by the current curl scaled by the gradient
                self.vx[i!(self.size,x,y)] += dt_s * self.vel_curl_at(x,y) * vort_grad_x;
                self.vy[i!(self.size,x,y)] += dt_s * self.vel_curl_at(x,y) * vort_grad_y;
            }
        }
    }

    fn vel_curl_at(&self, x: i32, y: i32) -> f32 {
        self.vx[i!(self.size, x,y+1)] - self.vx[i!(self.size, x,y-1)] +
        self.vy[i!(self.size, x-1,y)] - self.vy[i!(self.size, x+1,y)]
    }

    fn diffuse_field(&self, field: &[f32], rate: f32, is_vel: bool) -> Vec<f32> {
        let n = self.size;
        let mut diffused = vec![0.0; field.len()];
        for _ in 0..self.qual {
            let sol_0 = diffused.clone(); // This clone actually makes it faster?
            for y in 1..n-1 {
                for x in 1..n-1 {
                    // Diffusion is essentially a weighted averaging with neighbors
                    diffused[i!(n, x,y)] = 
                        (field[i!(n, x,y)] + rate * (
                            sol_0[i!(n, x,y+1)] + 
                            sol_0[i!(n, x,y-1)] + 
                            sol_0[i!(n, x+1,y)] + 
                            sol_0[i!(n, x-1,y)]
                        )) / (1.0 + 4.0 * rate);
                }
            }
            Fluid::bound_field(&mut diffused, n, is_vel);
        }
        diffused
    }

    fn advect_field(&self, field: &[f32], dt: f32) -> Vec<f32> {
        let n = self.size;
        let mut advected = vec![0.0; field.len()];
        for y in 1..n-1 {
            for x in 1..n-1 {
                // Trace velocities back to find value where it 'came' from
                let (vx, vy) = (self.vx[i!(n, x,y)], self.vy[i!(n, x,y)]);
                let (x0, y0) = (
                    (x as f32 - dt * vx).clamp(0.5, (n-1) as f32 - 0.5),
                    (y as f32 - dt * vy).clamp(0.5, (n-1) as f32 - 0.5));
                let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
                let s0 = field[i!(n, qx  , qy  )];
                let s1 = field[i!(n, qx+1, qy  )];
                let s2 = field[i!(n, qx  , qy+1)];
                let s3 = field[i!(n, qx+1, qy+1)];
                // Blerp between four samples to get proper value at position
                advected[i!(n, x ,y)] = blerp(s0, s1, s2, s3, x0.fract(), y0.fract());
            }
        }
        advected
    }

    fn bound_field(field: &mut [f32], n: i32, invert: bool) {
        let a = if invert { -1. } else { 1. };
        for i in 1..n-1 {
            field[i!(n,   0,i)] = a * field[i!(n,   1,i)];
            field[i!(n, n-1,i)] = a * field[i!(n, n-2,i)];
            field[i!(n, i,  0)] = a * field[i!(n, i,  1)];
            field[i!(n, i,n-1)] = a * field[i!(n, i,n-2)];
        }
    }
}