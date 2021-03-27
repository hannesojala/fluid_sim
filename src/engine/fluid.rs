// Window and simulation constants
pub const N: i32 = 256; 
pub const SIZE: usize = (N*N) as usize;
pub const RATIO: i32 = 2;
pub const MAX_FPS: u64 = 60;

// Quality related, 5 for low, 10-15 for medium, 20 for high
const GAUSS_SEIDEL_ITERATIONS: u32 = 10;

macro_rules! i(
    ($x:expr, $y:expr) => (
        (N * $y + $x) as usize
    )
);

fn bound_vfield(f: &mut Vec<(f32,f32)>) {
    for y in 1..N-1 {
        f[i!(0,y)].0 = -f[i!(1,y)].0;
        f[i!(N-1,y)].0 = -f[i!(N-2,y)].0;
    }
    for x in 1..N-1 {
        f[i!(x,0)].1 = -f[i!(x,1)].1;
        f[i!(x,N-1)].1 = -f[i!(x,N-2)].1;
    }
    f[i!(0, 0)].0 = (f[i!(0, 1)].0 + f[i!(1, 0)].0) / 2.0;
    f[i!(0, 0)].1 = (f[i!(0, 1)].1 + f[i!(1, 0)].1) / 2.0;

    f[i!(0, N-1)].0 = (f[i!(0, N-2)].0 + f[i!(1, N-1)].0) / 2.0;
    f[i!(0, N-1)].1 = (f[i!(0, N-2)].1 + f[i!(1, N-1)].1) / 2.0;

    f[i!(N-1, 0)].0 = (f[i!(N-1, 1)].0 + f[i!(N-2, 0)].0) / 2.0;
    f[i!(N-1, 0)].1 = (f[i!(N-1, 1)].1 + f[i!(N-2, 0)].1) / 2.0;

    f[i!(N-1, N-1)].0 = (f[i!(N-2, N-1)].0 + f[i!(N-1, N-2)].0) / 2.0;
    f[i!(N-1, N-1)].1 = (f[i!(N-2, N-1)].1 + f[i!(N-1, N-2)].1) / 2.0;
}

fn bound_sfield(f: &mut Vec<f32>) {
    for y in 1..N-1 {
        f[i!(0,y)] = -f[i!(1,y)];
        f[i!(N-1,y)] = -f[i!(N-2,y)];
    }
    for x in 1..N-1 {
        f[i!(x,0)] = -f[i!(x,1)];
        f[i!(x,N-1)] = -f[i!(x,N-2)];
    }
    f[i!(0, 0)] = (f[i!(0, 1)] + f[i!(1, 0)]) / 2.0;

    f[i!(0, N-1)] = (f[i!(0, N-2)] + f[i!(1, N-1)]) / 2.0;

    f[i!(N-1, 0)] = (f[i!(N-1, 1)] + f[i!(N-2, 0)]) / 2.0;

    f[i!(N-1, N-1)] = (f[i!(N-2, N-1)] + f[i!(N-1, N-2)]) / 2.0;
}

// Interpolations
fn lerp(v1: f32, v2: f32, k: f32) -> f32 {
    v1 + k * (v2 - v1)
}
// This is lazy
fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

// TODO: Generalize solve_field fns onto any dimension fields, not specific to diffusion
// Solves a scalar field as a linear system
fn solve_sfield(field: &Vec<f32>, a: f32) -> Vec<f32> {
    let f = field;
    let mut sol = vec![0.0; SIZE];
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        for y in 1..N-1 {
            for x in 1..N-1 {
                // Keep track of which cells are actually diffusable into
                // to prevent mass loss at walls. 
                let (mut n, mut s, mut e, mut w) 
                    = (1.,1.,1.,1.);
                    if x == 1 {  w = 0.0 }
                    if y == 1 {  n = 0.0 }
                    if x == N-2 {  e = 0.0 }
                    if y == N-2 {  s = 0.0 }
                sol[i!(x,y)] = 
                    (f[i!(x,y)] + a * (
                        s * sol[i!(x,y+1)] + 
                        n * sol[i!(x,y-1)] + 
                        e * sol[i!(x+1,y)] + 
                        w * sol[i!(x-1,y)]
                    )) / (1.0 + (n+s+e+w) * a);
            }
        }
        bound_sfield(&mut sol); // TODO: see how this works with my method above
    }
    sol
}
// Solves a vector field as a linear system
fn solve_vfield(field: &Vec<(f32,f32)>, a: f32) -> Vec<(f32,f32)> {
    let f = field;
    let mut sol = vec![(0.0, 0.0); SIZE]; // ewww
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        for y in 1..N-1 {
            for x in 1..N-1 {
                sol[i!(x,y)].0 = 
                    (f[i!(x,y)].0 + a * (
                        sol[i!(x,y+1)].0 + 
                        sol[i!(x,y-1)].0 + 
                        sol[i!(x+1,y)].0 + 
                        sol[i!(x-1,y)].0
                    )) / (1.0 + 4.0 * a);
                sol[i!(x,y)].1 = 
                    (f[i!(x,y)].1 + a * (
                        sol[i!(x,y+1)].1 + 
                        sol[i!(x,y-1)].1 + 
                        sol[i!(x+1,y)].1 + 
                        sol[i!(x-1,y)].1
                    )) / (1.0 + 4.0 * a);
            }
        }
        bound_vfield(&mut sol); // prevents diffusion into boundary
    }
    sol
}

// removes divergence from field, returns new field without divergence
// because helmholz theorem, all vector fields are a 
// sum of one with zero div and another with zero curl
fn enforce_div_eq_0(field: &Vec<(f32,f32)>) -> Vec<(f32,f32)> {
    let mut new = vec![(0.0, 0.0); SIZE];
    let mut div = vec![0.0; SIZE];
    let mut p = vec![0.0; SIZE];
    let h = 1.0 / N as f32;
    for y in 1..N-1 {
        for x in 1..N-1 {
            div[i!(x,y)] = -0.5 * h * (field[i!(x+1,y)].0 - field[i!(x-1,y)].0 + field[i!(x,y+1)].1 - field[i!(x,y-1)].1);
        }
    }
    bound_sfield(&mut div);
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        for y in 1..N-1 {
            for x in 1..N-1 {
                p[i!(x,y)] = (div[i!(x,y)]+p[i!(x-1,y)]+p[i!(x+1,y)]+p[i!(x,y-1)]+p[i!(x,y+1)]) / 4.0;
            }
        }
        bound_sfield(&mut p);
    }
    for y in 1..N-1 {
        for x in 1..N-1 {
            new[i!(x,y)].0 = field[i!(x,y)].0 - (0.5*(p[i!(x+1,y)] - p[i!(x-1,y)]) / h);
            new[i!(x,y)].1 = field[i!(x,y)].1 - (0.5*(p[i!(x,y+1)] - p[i!(x,y-1)]) / h);
        }
    }
    bound_vfield(&mut new);
    new
}

// advects vector field along itself
fn advect_vfield(field: &Vec<(f32,f32)>, dt: f32) -> Vec<(f32,f32)> {
    let mut new = vec![(0.0,0.0); SIZE];
    for y in 1..N-1 {
        for x in 1..N-1 {
            let (vx, vy) = field[i!(x,y)];
            let (x0, y0) = (
                (x as f32 - dt * vx).clamp(0.5, (N-1) as f32 - 0.5),
                (y as f32 - dt * vy).clamp(0.5, (N-1) as f32 - 0.5));
            let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
            let (v0x, v0y) = field[i!(qx  , qy  )];
            let (v1x, v1y) = field[i!(qx+1, qy  )];
            let (v2x, v2y) = field[i!(qx  , qy+1)];
            let (v3x, v3y) = field[i!(qx+1, qy+1)];
            new[i!(x ,y)] = (blerp(v0x, v1x, v2x, v3x, x0.fract(), y0.fract()),
                             blerp(v0y, v1y, v2y, v3y, x0.fract(), y0.fract()));
        }
    }
    bound_vfield(&mut new);
    new
}

// advects a scalar field along a vector field
fn advect_sfield(sfield: &Vec<f32>, vfield: &Vec<(f32,f32)>, dt: f32) -> Vec<f32> {
    let mut new = vec![0.0; SIZE];
    for y in 1..N-1 {
        for x in 1..N-1 {
            let (vx, vy) = vfield[i!(x,y)];
            let (x0, y0) = (
                (x as f32 - dt * vx).clamp(0.5, (N-1) as f32 - 0.5),
                (y as f32 - dt * vy).clamp(0.5, (N-1) as f32 - 0.5));
            let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
            let s0 = sfield[i!(qx  , qy  )];
            let s1 = sfield[i!(qx+1, qy  )];
            let s2 = sfield[i!(qx  , qy+1)];
            let s3 = sfield[i!(qx+1, qy+1)];
            new[i!(x ,y)] = blerp(s0, s1, s2, s3, x0.fract(), y0.fract());
        }
    }
    bound_sfield(&mut new);
    new
}

pub struct Fluid {
    pub visc: f32, pub diff: f32,
    vel: Vec<(f32,f32) >,
    dye: Vec<f32>
}

impl Fluid {
    pub fn new(visc: f32, diff: f32) -> Fluid {
        Fluid {
            visc, diff,
            vel:  vec![(0.0,0.0); SIZE],
            dye:  vec![0.0; SIZE],
        }
    }

    pub fn update(&mut self, dt_s: f32) {
        // let mass_before: f32 = self.dye.iter().sum(); // see below
        self.vel = solve_vfield(&self.vel, dt_s * self.visc * ((N - 2)*(N - 2)) as f32);
        self.vel =  enforce_div_eq_0(&self.vel);

        self.vel = advect_vfield(&self.vel, dt_s);
        self.vel =  enforce_div_eq_0(&self.vel);

        self.dye = solve_sfield(&self.dye, dt_s * self.diff * ((N - 2)*(N - 2)) as f32);
        self.dye = advect_sfield(&self.dye, &self.vel, dt_s);
        // let mass_after: f32 = self.dye.iter().sum();
        // fails sometimes due to advection
        // works with diffusion
        // todo: enforce mass conservation of dye during advection
        // assert!((mass_before - mass_after).abs() <= 0.01 * mass_before); // Make sure dye mass is conserved
    }

    pub fn add_dye(&mut self, x: i32, y: i32, amt: f32) {
        self.dye[i!(x, y)] += amt;
    }

    pub fn set_dye(&mut self, x: i32, y: i32, amt: f32) {
        self.dye[i!(x, y)] = amt;
    }

    pub fn add_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        self.vel[i!(x, y)].0 += vx;
        self.vel[i!(x, y)].1 += vy;
    }

    pub fn set_vel(&mut self, x: i32, y: i32, vx: f32, vy: f32) {
        self.vel[i!(x, y)].0 = vx;
        self.vel[i!(x, y)].1 = vy;
    }

    pub fn dye(&mut self, x: i32, y: i32) -> f32 {
        self.dye[i!(x,y)]
    }

    pub fn vel(&mut self, x: i32, y: i32) -> (f32,f32) {
        self.vel[i!(x,y)]
    }
}