// Window and simulation constants
pub const N: i32 = 124; 
pub const SIZE: usize = (N*N) as usize;
pub const RATIO: i32 = 4;
pub const MAX_FPS: u64 = 300;

pub const DIFF: f32 = 1e-4; // What should these be???
pub const VISC: f32 = 1e-4;

const GAUSS_SEIDEL_ITERATIONS: u32 = 4;

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

pub fn solve_sfield(field: &Vec<f32>, a: f32) -> Vec<f32> {
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
        
    }
    sol
}

// Solves a 2d vector field as a linear system with coeffs a and b
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
        bound_vel(&mut sol); // prevents diffusion into boundary
    }
    sol
}



pub fn bound_vel(f: &mut Vec<(f32,f32)>) {
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