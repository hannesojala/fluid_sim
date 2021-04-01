// Quality related, 5 for low, 10-15 for medium, 20 for high
const GAUSS_SEIDEL_ITERATIONS: u32 = 15;

// TODO: Use rayon to parallelize?

macro_rules! i(
    ($n:expr, $x:expr, $y:expr) => (
        ($n * $y + $x) as usize
    )
);

// boundary conditions
fn bound(field: &mut Vec<f32>, n: i32, invert: bool) {
    let a = if invert { -1. } else { 1. };
    for i in 1..n-1 {
        field[i!(n,   0,i)] = a * field[i!(n,   1,i)];
        field[i!(n, n-1,i)] = a * field[i!(n, n-2,i)];
        field[i!(n, i,  0)] = a * field[i!(n, i,  1)];
        field[i!(n, i,n-1)] = a * field[i!(n, i,n-2)];
    }
}

// Interpolations
fn lerp(v1: f32, v2: f32, k: f32) -> f32 { v1 + k * (v2 - v1) }

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

// diffuses a scalar field at a rate. is_vel is for boundary conditions
fn diffuse(field: &Vec<f32>, a: f32, n: i32, is_vel: bool) -> Vec<f32> {
    let mut sol = vec![0.0; field.len()];
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        let sol_0 = sol.clone(); // This clone actually makes it faster?
        for y in 1..n-1 {
            for x in 1..n-1 {
                sol[i!(n, x,y)] = 
                    (field[i!(n, x,y)] + a * (
                        sol_0[i!(n, x,y+1)] + 
                        sol_0[i!(n, x,y-1)] + 
                        sol_0[i!(n, x+1,y)] + 
                        sol_0[i!(n, x-1,y)]
                    )) / (1.0 + 4.0 * a);
            }
        }
        bound(&mut sol, n, is_vel);
    }
    sol
}

// removes divergence
fn remove_div(fieldx: &mut Vec<f32>, fieldy: &mut Vec<f32>, n: i32) {
    let mut div = vec![0.0; fieldx.len()];
    let mut p = vec![0.0; fieldx.len()];
    for y in 1..n-1 {
        for x in (1..n-1).into_iter() {
            div[i!(n, x,y)] = -0.5 * (fieldx[i!(n, x+1,y)] - fieldx[i!(n, x-1,y)] + fieldy[i!(n, x,y+1)] - fieldy[i!(n, x,y-1)]) / n as f32;
        }
    }
    bound(&mut div, n, false);
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        let p_0 = p.clone(); // This clone actually makes it faster?
        for y in 1..n-1 {
            for x in 1..n-1 {
                p[i!(n, x,y)] = (div[i!(n, x,y)] + p_0[i!(n, x-1,y)] + p_0[i!(n, x+1,y)] + p_0[i!(n, x,y-1)] + p_0[i!(n, x,y+1)]) / 4.0;
            }
        }
        bound(&mut p, n, false);
    }
    for y in 1..n-1 {
        for x in 1..n-1 {
            fieldx[i!(n, x,y)] -= 0.5 * (p[i!(n, x+1,y)] - p[i!(n, x-1,y)]) * n as f32;
            fieldy[i!(n, x,y)] -= 0.5 * (p[i!(n, x,y+1)] - p[i!(n, x,y-1)]) * n as f32;
        }
    }
    bound(fieldx, n, true);
    bound(fieldy, n, true);
}

// advects a scalar field along the velocity field
fn advect(sfield: & Vec<f32>, fluid: & Fluid, dt: f32, is_vel: bool) -> Vec<f32> {
    let n = fluid.size;
    let mut new = vec![0.0; sfield.len()];
    for y in 1..n-1 {
        for x in 1..n-1 {
            let (vx, vy) = (fluid.vx[i!(n, x,y)], fluid.vy[i!(n, x,y)]);
            let (x0, y0) = (
                (x as f32 - dt * vx).clamp(0.5, (n-1) as f32 - 0.5),
                (y as f32 - dt * vy).clamp(0.5, (n-1) as f32 - 0.5));
            let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
            let s0 = sfield[i!(n, qx  , qy  )];
            let s1 = sfield[i!(n, qx+1, qy  )];
            let s2 = sfield[i!(n, qx  , qy+1)];
            let s3 = sfield[i!(n, qx+1, qy+1)];
            new[i!(n, x ,y)] = blerp(s0, s1, s2, s3, x0.fract(), y0.fract());
        }
    }
    bound(&mut new, n, is_vel);
    new
}

pub struct Fluid {
    pub visc: f32,
    pub diff: f32,
    size: i32,
    vx: Vec<f32>,
    vy: Vec<f32>,
    dye: Vec<f32>
}

impl Fluid {
    pub fn new(visc: f32, diff: f32, size: i32) -> Fluid {
        let flat_size = ((size+2)*(size+2)) as usize;
        Fluid {
            visc, diff, size,
            vx:  vec![0.0; flat_size],
            vy:  vec![0.0; flat_size],
            dye:  vec![0.0; flat_size],
        }
    }

    pub fn update(&mut self, dt_s: f32) {
        self.vx = diffuse(&mut self.vx, dt_s * self.visc * ((self.size - 2)*(self.size - 2)) as f32, self.size, true);
        self.vy = diffuse(&mut self.vy, dt_s * self.visc * ((self.size - 2)*(self.size - 2)) as f32, self.size, true);
        
        remove_div(&mut self.vx, &mut self.vy, self.size);

        self.vx = advect(& self.vx, & self, dt_s, true);
        self.vy = advect(& self.vy, & self, dt_s, true);
        remove_div(&mut self.vx, &mut self.vy, self.size);
        
        diffuse(&mut self.dye, dt_s * self.diff * ((self.size - 2)*(self.size - 2)) as f32, self.size, false);
        self.dye = advect(& self.dye, &self, dt_s, false);
    }

    pub fn add_dye(&mut self, x: i32, y: i32, amt: f32) {
        self.dye[i!(self.size, x, y)] += amt;
    }

    pub fn add_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        self.vx[i!(self.size, x, y)] += vx;
        self.vy[i!(self.size, x, y)] += vy;
    }

    pub fn dye_at(&mut self, x: i32, y: i32) -> f32 {
        self.dye[i!(self.size, x,y)]
    }

    pub fn vel_at(&mut self, x: i32, y: i32) -> (f32,f32) {
        (self.vx[i!(self.size, x,y)], self.vy[i!(self.size, x,y)])
    }
}