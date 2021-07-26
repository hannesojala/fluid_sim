// Hate this code so much but IDK how to make it prettier.

// Quality related, 5 for low, 10-15 for medium, 20 for high
const GAUSS_SEIDEL_ITERATIONS: u32 = 20;

// TODO: Use rayon to parallelize?

macro_rules! i(
    ($n:expr, $x:expr, $y:expr) => (
        ($n * $y + $x) as usize
    )
);

fn lerp(v1: f32, v2: f32, k: f32) -> f32 { v1 + k * (v2 - v1) }

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

pub struct Fluid {
    pub visc: f32,
    pub diff: f32,
    size: i32,
    vx: Vec<f32>,
    vy: Vec<f32>,
    dye_r: Vec<f32>,
    dye_g: Vec<f32>,
    dye_b: Vec<f32>,
}

impl Fluid {
    pub fn new(visc: f32, diff: f32, size: i32) -> Fluid {
        let flat_size = ((size+2)*(size+2)) as usize;
        Fluid {
            visc, diff, size,
            vx:  vec![0.0; flat_size],
            vy:  vec![0.0; flat_size],
            dye_r:  vec![0.0; flat_size],
            dye_g:  vec![0.0; flat_size],
            dye_b:  vec![0.0; flat_size],
        }
    }

    pub fn update(&mut self, dt_s: f32) {
        let velocity_diffusion_rate = dt_s * self.visc * ((self.size - 2)*(self.size - 2)) as f32;
        self.vx = Fluid::diffuse(&self.vx, velocity_diffusion_rate, self.size, true);
        self.vy = Fluid::diffuse(&self.vy, velocity_diffusion_rate, self.size, true);
        self.remove_div();

        self.vx = Fluid::advect(&self.vx, &self.vx, &self.vy, dt_s, self.size);
        Fluid::bound(&mut self.vx, self.size, true);
        self.vy = Fluid::advect(&self.vy, &self.vx, &self.vy, dt_s, self.size);
        Fluid::bound(&mut self.vx, self.size, true);
        self.remove_div();
        
        let dye_diffusion_rate = dt_s * self.diff * ((self.size - 2)*(self.size - 2)) as f32;
        self.dye_r = Fluid::diffuse(&self.dye_r, dye_diffusion_rate, self.size, false);
        self.dye_r = Fluid::advect(&self.dye_r, &self.vx, &self.vy, dt_s, self.size);
        Fluid::bound(&mut self.dye_r, self.size, false);

        self.dye_g = Fluid::diffuse(&self.dye_g, dye_diffusion_rate, self.size, false);
        self.dye_g = Fluid::advect(&self.dye_g, &self.vx, &self.vy, dt_s, self.size);
        Fluid::bound(&mut self.dye_g, self.size, false);

        self.dye_b = Fluid::diffuse(&self.dye_b, dye_diffusion_rate, self.size, false);
        self.dye_b = Fluid::advect(&self.dye_b, &self.vx, &self.vy, dt_s, self.size);
        Fluid::bound(&mut self.dye_b, self.size, false);
    }

    pub fn add_dye(&mut self, x: i32, y: i32, rgb: (f32, f32, f32)) {
        self.dye_r[i!(self.size, x, y)] += rgb.0;
        self.dye_g[i!(self.size, x, y)] += rgb.1;
        self.dye_b[i!(self.size, x, y)] += rgb.2;
    }

    pub fn add_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        self.vx[i!(self.size, x, y)] += vx;
        self.vy[i!(self.size, x, y)] += vy;
    }

    pub fn dye_at(&mut self, x: i32, y: i32) -> (f32, f32, f32) {
        let i = i!(self.size, x,y);
        (self.dye_r[i], self.dye_g[i], self.dye_b[i])
    }

    pub fn vel_at(&mut self, x: i32, y: i32) -> (f32,f32) {
        (self.vx[i!(self.size, x,y)], self.vy[i!(self.size, x,y)])
    }

    fn remove_div(&mut self) {
        let n = self.size;
        let mut div = vec![0.0; self.vx.len()];
        let mut p = vec![0.0; self.vx.len()];
        for y in 1..n-1 {
            for x in 1..n-1 {
                div[i!(n, x,y)] = -0.5 * (self.vx[i!(n, x+1,y)] - self.vx[i!(n, x-1,y)] + self.vy[i!(n, x,y+1)] - self.vy[i!(n, x,y-1)]) / n as f32;
            }
        }
        Fluid::bound(&mut div, n, false);
        for _ in 0..GAUSS_SEIDEL_ITERATIONS {
            let p_0 = p.clone(); // This clone actually makes it faster?
            for y in 1..n-1 {
                for x in 1..n-1 {
                    p[i!(n, x,y)] = (div[i!(n, x,y)] + p_0[i!(n, x-1,y)] + p_0[i!(n, x+1,y)] + p_0[i!(n, x,y-1)] + p_0[i!(n, x,y+1)]) / 4.0;
                }
            }
            Fluid::bound(&mut p, n, false);
        }
        for y in 1..n-1 {
            for x in 1..n-1 {
                self.vx[i!(n, x,y)] -= 0.5 * (p[i!(n, x+1,y)] - p[i!(n, x-1,y)]) * n as f32;
                self.vy[i!(n, x,y)] -= 0.5 * (p[i!(n, x,y+1)] - p[i!(n, x,y-1)]) * n as f32;
            }
        }
        Fluid::bound(&mut self.vx, n, true);
        Fluid::bound(&mut self.vy, n, true);
    }

    fn diffuse(field: &[f32], rate: f32, n: i32, is_vel: bool) -> Vec<f32> {
        let mut sol = vec![0.0; field.len()];
        for _ in 0..GAUSS_SEIDEL_ITERATIONS {
            let sol_0 = sol.clone(); // This clone actually makes it faster?
            for y in 1..n-1 {
                for x in 1..n-1 {
                    sol[i!(n, x,y)] = 
                        (field[i!(n, x,y)] + rate * (
                            sol_0[i!(n, x,y+1)] + 
                            sol_0[i!(n, x,y-1)] + 
                            sol_0[i!(n, x+1,y)] + 
                            sol_0[i!(n, x-1,y)]
                        )) / (1.0 + 4.0 * rate);
                }
            }
            Fluid::bound(&mut sol, n, is_vel);
        }
        sol
    }

    fn advect(field: &[f32], vx: &[f32], vy: &[f32], dt: f32, n: i32) -> Vec<f32> {
        let mut new = vec![0.0; field.len()];
        for y in 1..n-1 {
            for x in 1..n-1 {
                let (vx, vy) = (vx[i!(n, x,y)], vy[i!(n, x,y)]);
                let (x0, y0) = (
                    (x as f32 - dt * vx).clamp(0.5, (n-1) as f32 - 0.5),
                    (y as f32 - dt * vy).clamp(0.5, (n-1) as f32 - 0.5));
                let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
                let s0 = field[i!(n, qx  , qy  )];
                let s1 = field[i!(n, qx+1, qy  )];
                let s2 = field[i!(n, qx  , qy+1)];
                let s3 = field[i!(n, qx+1, qy+1)];
                new[i!(n, x ,y)] = blerp(s0, s1, s2, s3, x0.fract(), y0.fract());
            }
        }
        new
    }

    fn bound(field: &mut [f32], n: i32, invert: bool) {
        let a = if invert { -1. } else { 1. };
        for i in 1..n-1 {
            field[i!(n,   0,i)] = a * field[i!(n,   1,i)];
            field[i!(n, n-1,i)] = a * field[i!(n, n-2,i)];
            field[i!(n, i,  0)] = a * field[i!(n, i,  1)];
            field[i!(n, i,n-1)] = a * field[i!(n, i,n-2)];
        }
    }
}