use std::thread;

// Quality related, 5 for low, 10-15 for medium, 20 for high
const GAUSS_SEIDEL_ITERATIONS: u32 = 20;
// Main performance bottleneck is the solving of linear eqs 
// for diffusion and divergence.

// TODO: rewrite with iterators and use rayon to parallelize?

macro_rules! i(
    ($n:expr, $x:expr, $y:expr) => (
        ($n * $y + $x) as usize
    )
);

fn bound_vel(field: &mut Vec<(f32,f32)>, N: i32) {
    for i in 1..N-1 {
        field[i!(N,   0,i)].0 = -field[i!(N, 1,  i)].0;
        field[i!(N, N-1,i)].0 = -field[i!(N, N-2,i)].0;
        field[i!(N, i,  0)].1 = -field[i!(N, i,  1)].1;
        field[i!(N, i,N-1)].1 = -field[i!(N, i,N-2)].1;
    }
}

fn bound_scalar(field: &mut Vec<f32>, N: i32) {
    for i in 1..N-1 {
        field[i!(N,   0,i)] = field[i!(N,   1,i)];
        field[i!(N, N-1,i)] = field[i!(N, N-2,i)];
        field[i!(N, i,  0)] = field[i!(N, i,  1)];
        field[i!(N, i,N-1)] = field[i!(N, i,N-2)];
    }
}

// Interpolations
fn lerp(v1: f32, v2: f32, k: f32) -> f32 { v1 + k * (v2 - v1) }

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

// TODO: Generalize solve_field fns onto any dimension fields, not specific to diffusion
fn diffuse_dye(field: &mut Vec<f32>, a: f32, N: i32) {
    let mut sol = vec![0.0; field.len()];
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        let sol_0 = sol.clone(); // This clone actually makes it faster?
        for y in 1..N-1 {
            for x in 1..N-1 {
                sol[i!(N, x,y)] = 
                    (field[i!(N, x,y)] + a * (
                        sol_0[i!(N, x,y+1)] + 
                        sol_0[i!(N, x,y-1)] + 
                        sol_0[i!(N, x+1,y)] + 
                        sol_0[i!(N, x-1,y)]
                    )) / (1.0 + 4.0 * a);
            }
        }
        bound_scalar(&mut sol, N);
    }
    *field = sol;
}

fn diffuse_vel(field: &mut Vec<(f32,f32)>, a: f32, N: i32) {
    let mut sol = vec![(0.0, 0.0); field.len()];
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        let sol_0 = sol.clone(); // This clone actually makes it faster?
        for y in 1..N-1 {
            for x in 1..N-1 {
                sol[i!(N, x,y)].0 = 
                    (field[i!(N, x,y)].0 + a * (
                        sol_0[i!(N, x,y+1)].0 + 
                        sol_0[i!(N, x,y-1)].0 + 
                        sol_0[i!(N, x+1,y)].0 + 
                        sol_0[i!(N, x-1,y)].0
                    )) / (1.0 + 4.0 * a);
                sol[i!(N, x,y)].1 = 
                    (field[i!(N, x,y)].1 + a * (
                        sol_0[i!(N, x,y+1)].1 + 
                        sol_0[i!(N, x,y-1)].1 + 
                        sol_0[i!(N, x+1,y)].1 + 
                        sol_0[i!(N, x-1,y)].1
                    )) / (1.0 + 4.0 * a);
            }
        }
        bound_vel(&mut sol, N);
    }
    *field = sol;
}

fn remove_div(field: &mut Vec<(f32,f32)>, N: i32) {
    let mut div = vec![0.0; field.len()];
    let mut p = vec![0.0; field.len()];
    for y in 1..N-1 {
        for x in 1..N-1 {
            div[i!(N, x,y)] = -0.5 * (field[i!(N, x+1,y)].0 - field[i!(N, x-1,y)].0 + field[i!(N, x,y+1)].1 - field[i!(N, x,y-1)].1) / N as f32;
        }
    }
    bound_scalar(&mut div, N);
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        let p_0 = p.clone(); // This clone actually makes it faster?
        for y in 1..N-1 {
            for x in 1..N-1 {
                p[i!(N, x,y)] = (div[i!(N, x,y)] + p_0[i!(N, x-1,y)] + p_0[i!(N, x+1,y)] + p_0[i!(N, x,y-1)] + p_0[i!(N, x,y+1)]) / 4.0;
            }
        }
        bound_scalar(&mut p, N);
    }
    for y in 1..N-1 {
        for x in 1..N-1 {
            field[i!(N, x,y)].0 -= 0.5 * (p[i!(N, x+1,y)] - p[i!(N, x-1,y)]) * N as f32;
            field[i!(N, x,y)].1 -= 0.5 * (p[i!(N, x,y+1)] - p[i!(N, x,y-1)]) * N as f32;
        }
    }
    bound_vel(field, N);
}

// advects vector field along itself
fn advect_vel(field: &mut Vec<(f32,f32)>, dt: f32, N: i32) {
    let mut new = vec![(0.0,0.0); field.len()];
    for y in 1..N-1 {
        for x in 1..N-1 {
            let (vx, vy) = field[i!(N, x,y)];
            let (x0, y0) = (
                (x as f32 - dt * vx).clamp(0.5, (N-1) as f32 - 0.5),
                (y as f32 - dt * vy).clamp(0.5, (N-1) as f32 - 0.5));
            let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
            let (v0x, v0y) = field[i!(N, qx  , qy  )];
            let (v1x, v1y) = field[i!(N, qx+1, qy  )];
            let (v2x, v2y) = field[i!(N, qx  , qy+1)];
            let (v3x, v3y) = field[i!(N, qx+1, qy+1)];
            new[i!(N, x ,y)] = (blerp(v0x, v1x, v2x, v3x, x0.fract(), y0.fract()),
                             blerp(v0y, v1y, v2y, v3y, x0.fract(), y0.fract()));
        }
    }
    bound_vel(&mut new, N);
    *field = new;
}

// advects a scalar field along a vector field
fn advect_dye(sfield: &mut Vec<f32>, vfield: &Vec<(f32,f32)>, dt: f32, N: i32) {
    let mut new = vec![0.0; sfield.len()];
    for y in 1..N-1 {
        for x in 1..N-1 {
            let (vx, vy) = vfield[i!(N, x,y)];
            let (x0, y0) = (
                (x as f32 - dt * vx).clamp(0.5, (N-1) as f32 - 0.5),
                (y as f32 - dt * vy).clamp(0.5, (N-1) as f32 - 0.5));
            let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
            let s0 = sfield[i!(N, qx  , qy  )];
            let s1 = sfield[i!(N, qx+1, qy  )];
            let s2 = sfield[i!(N, qx  , qy+1)];
            let s3 = sfield[i!(N, qx+1, qy+1)];
            new[i!(N, x ,y)] = blerp(s0, s1, s2, s3, x0.fract(), y0.fract());
        }
    }
    bound_scalar(&mut new, N);
    *sfield = new;
}

pub struct Fluid {
    pub visc: f32,
    pub diff: f32,
    size: i32,
    vel: Vec<(f32,f32) >,
    dye: Vec<f32>
}

impl Fluid {
    pub fn new(visc: f32, diff: f32, size: i32) -> Fluid {
        let flat_size = ((size+2)*(size+2)) as usize;
        Fluid {
            visc, diff, size,
            vel:  vec![(0.0,0.0); flat_size],
            dye:  vec![0.0; flat_size],
        }
    }

    pub fn update(&mut self, dt_s: f32) {
        diffuse_vel(&mut self.vel, dt_s * self.visc * ((self.size - 2)*(self.size - 2)) as f32, self.size);
        remove_div(&mut self.vel, self.size);

        advect_vel(&mut self.vel, dt_s, self.size);
        remove_div(&mut self.vel, self.size);

        diffuse_dye(&mut self.dye, dt_s * self.diff * ((self.size - 2)*(self.size - 2)) as f32, self.size);
        advect_dye(&mut self.dye, &self.vel, dt_s, self.size);
    }

    pub fn add_dye(&mut self, x: i32, y: i32, amt: f32) {
        self.dye[i!(self.size, x, y)] += amt;
    }

    pub fn add_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        self.vel[i!(self.size, x, y)].0 += vx;
        self.vel[i!(self.size, x, y)].1 += vy;
    }

    pub fn dye(&mut self, x: i32, y: i32) -> f32 {
        self.dye[i!(self.size, x,y)]
    }

    pub fn vel(&mut self, x: i32, y: i32) -> (f32,f32) {
        self.vel[i!(self.size, x,y)]
    }
}