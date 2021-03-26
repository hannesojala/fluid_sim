// Window and simulation constants
pub const N: i32 = 256; 
pub const SIZE: usize = (N*N) as usize;
pub const RATIO: i32 = 2;
pub const MAX_FPS: u64 = 300;

pub const DIFF: f32 = 1e-6; // What should these be???
pub const VISC: f32 = 1e-6;

const GAUSS_SEIDEL_ITERATIONS: u32 = 10;

macro_rules! i(
    ($x:expr, $y:expr) => (
        (N * $y + $x) as usize
    )
);

// Interpolations
pub fn lerp(v1: f32, v2: f32, k: f32) -> f32 {
    v1 + k * (v2 - v1)
}

pub fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

// TODO: Generalize solve_field fns onto any dimension fields
// Solves a scalar field as a linear system with coeffs a and b
pub fn solve_sfield(field: &Vec<f32>, a: f32) -> Vec<f32> {
    let f = field;
    let mut sol = vec![0.0; SIZE];
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        for y in 1..N-1 {
            for x in 1..N-1 {
                // Keep track of which cells are actually diffusable into
                // to prevent mass loss at walls. 
                // does bound_sfield replace this?
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
        
    }
    sol
}

pub fn solve_vfield(field: &Vec<(f32,f32)>, a: f32) -> Vec<(f32,f32)> {
    let f = field;
    let mut sol = vec![(0.0, 0.0); SIZE]; // ewww
    for _ in 0..GAUSS_SEIDEL_ITERATIONS {
        for y in 1..N-1 {
            for x in 1..N-1 {
                // Should I prevent velocity loss at walls?
                // Other implementations set the boundary velocity to the edge velocity
                // with reversed normal component
                // Only that component erased, maybe?
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

pub fn enforce_div_eq_0(field: &Vec<(f32,f32)>) -> Vec<(f32,f32)> {
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
pub fn advect_vfield(field: &Vec<(f32,f32)>, dt: f32) -> Vec<(f32,f32)> {
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
pub fn advect_sfield(sfield: &Vec<f32>, vfield: &Vec<(f32,f32)>, dt: f32) -> Vec<f32> {
    let mut new = vec![0.0; SIZE];
    for y in 1..N-1 {
        for x in 1..N-1 {
            let s = sfield[i!(x,y)];
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

pub fn bound_vfield(f: &mut Vec<(f32,f32)>) {
    // Walls, unchecked
    for y in 1..N-1 {
        f[i!(0,y)].0 = -f[i!(1,y)].0;
        f[i!(N-1,y)].0 = -f[i!(N-2,y)].0;
    }
    for x in 1..N-1 {
        f[i!(x,0)].1 = -f[i!(x,1)].1;
        f[i!(x,N-1)].1 = -f[i!(x,N-2)].1;
    }
    // Corners, unchecked
    f[i!(0, 0)].0 = (f[i!(0, 1)].0 + f[i!(1, 0)].0) / 2.0;
    f[i!(0, 0)].1 = (f[i!(0, 1)].1 + f[i!(1, 0)].1) / 2.0;

    f[i!(0, N-1)].0 = (f[i!(0, N-2)].0 + f[i!(1, N-1)].0) / 2.0;
    f[i!(0, N-1)].1 = (f[i!(0, N-2)].1 + f[i!(1, N-1)].1) / 2.0;

    f[i!(N-1, 0)].0 = (f[i!(N-1, 1)].0 + f[i!(N-2, 0)].0) / 2.0;
    f[i!(N-1, 0)].1 = (f[i!(N-1, 1)].1 + f[i!(N-2, 0)].1) / 2.0;

    f[i!(N-1, N-1)].0 = (f[i!(N-2, N-1)].0 + f[i!(N-1, N-2)].0) / 2.0;
    f[i!(N-1, N-1)].1 = (f[i!(N-2, N-1)].1 + f[i!(N-1, N-2)].1) / 2.0;
}

pub fn bound_sfield(f: &mut Vec<f32>) {
    // Walls, unchecked
    for y in 1..N-1 {
        f[i!(0,y)] = -f[i!(1,y)];
        f[i!(N-1,y)] = -f[i!(N-2,y)];
    }
    for x in 1..N-1 {
        f[i!(x,0)] = -f[i!(x,1)];
        f[i!(x,N-1)] = -f[i!(x,N-2)];
    }
    // Corners, unchecked
    f[i!(0, 0)] = (f[i!(0, 1)] + f[i!(1, 0)]) / 2.0;

    f[i!(0, N-1)] = (f[i!(0, N-2)] + f[i!(1, N-1)]) / 2.0;

    f[i!(N-1, 0)] = (f[i!(N-1, 1)] + f[i!(N-2, 0)]) / 2.0;

    f[i!(N-1, N-1)] = (f[i!(N-2, N-1)] + f[i!(N-1, N-2)]) / 2.0;
}