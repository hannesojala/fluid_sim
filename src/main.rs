use std::thread::sleep;
use std::time::{Duration, Instant};
use sdl2::event::Event;
use sdl2::pixels::Color;
use sdl2::rect::Rect;

// TODO: Fix code, simplify, figure out vorticity
// Clean shit
// learn about macros
// WHERE IS MA VORTICITEEE

const SIZE: u32 = 128;
const RATIO: u32 = 4;
const WINDOW_SIZE: u32 = SIZE * RATIO;
const FLAT_SIZE: usize = (SIZE * SIZE) as usize;
const MAX_FPS: u64 = 144;
const MIN_DT : Duration = Duration::from_millis(1000 / MAX_FPS);

const DYE_DIFF_RATE: f32 = 0.005;
const VISC_THINGY: f32 = 1.0;
const GAUSS_SEIDEL_STEPS: u32 = 8;

fn initial_velocity_fn(x:f32, y:f32) -> (f32, f32) {
    (0.0,0.0)//(8. * -y, 8. * x)
}

fn coord_to_idx((x, y): (i32,i32)) -> usize {
    (SIZE as i32 * y.rem_euclid(SIZE as i32) + x.rem_euclid(SIZE as i32)) as usize
}

fn idx_to_coord(i: usize) -> (u32, u32) {
    (i as u32 % SIZE, i as u32 / SIZE)
}

fn lerp(v1: f32, v2: f32, k: f32) -> f32 {
    v1 + k * (v2 - v1)
}

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 {
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
}

fn field_to_canvas_center((x, y): (u32, u32)) -> (i32, i32) {
    let (x, y) = (RATIO * x + RATIO / 2, RATIO * y + RATIO / 2);
    (x as i32, y as i32)
}

fn field_to_canvas_corner((x, y): (u32, u32)) -> (i32, i32) {
    let (x, y) = (RATIO * x, RATIO * y);
    (x as i32, y as i32)
}

fn main() {
    let mut engine = Engine::init();
    while engine.running {
        engine.update();
        engine.events();
        engine.render();
    }
}

struct Engine {
    canvas: sdl2::render::WindowCanvas,
    event_pump: sdl2::EventPump,
    running: bool,
    current_time: Instant,
    dye_field: Vec<f32>,
    velocity_field: Vec<(f32, f32)>
}

impl Engine {
    fn init() -> Engine {
        let sdl_context = sdl2::init().unwrap();
        let video_subsystem = sdl_context.video().unwrap();
        let wsize = WINDOW_SIZE;
        let window = video_subsystem.window("fluid", wsize, wsize)
            .position_centered()
            .build().unwrap();

        let dye_field = vec![0.0; FLAT_SIZE];

        let mut velocity_field = vec![(0.0,0.0); FLAT_SIZE];
        for (i, vel) in velocity_field.iter_mut().enumerate() {
            let pt = idx_to_coord(i);
            let (x, y) = (pt.0 as f32 - SIZE as f32 / 2.0, pt.1 as f32 - SIZE as f32 / 2.0);
            *vel = initial_velocity_fn(x, y);
        }
        
        Engine {
            canvas: window.into_canvas().build().unwrap(),
            event_pump: sdl_context.event_pump().unwrap(),
            running: true,
            current_time: Instant::now(),
            dye_field,
            velocity_field
        }
    }

    fn diffuse_dye(&mut self, dt: f32) { // may be wrong
        let rate = dt * DYE_DIFF_RATE;// * (SIZE as f32 * SIZE as f32);
        for k in 0..GAUSS_SEIDEL_STEPS { // is 20 in paper
            let mut updated_dye_field = vec![0.0; FLAT_SIZE];
            for (i, dye) in self.dye_field.iter().enumerate() {
                let (x, y) = idx_to_coord(i);
                updated_dye_field[i] = (self.dye_field[i] + rate * (
                    self.dye_field[coord_to_idx((x as i32   , y as i32 +1))] +
                    self.dye_field[coord_to_idx((x as i32   , y as i32 -1))] +
                    self.dye_field[coord_to_idx((x as i32 +1, y as i32   ))] +
                    self.dye_field[coord_to_idx((x as i32 -1, y as i32   ))]
                )) / (1.0 + 4.0 * rate);
            }
            self.dye_field = updated_dye_field;
        }
    }

    fn advect_dye(&mut self, dt: f32) {
        let mut updated_dye_field = vec![0.0; FLAT_SIZE];
        for (i, (vx, vy)) in self.velocity_field.iter().enumerate() {
            let (x, y) = idx_to_coord(i);
            let (x_pre, y_pre) = (
                (x as f32 - dt * vx).rem_euclid(SIZE as f32),
                (y as f32 - dt * vy).rem_euclid(SIZE as f32));
            let (qx, qy) = (x_pre.floor() as i32 , y_pre.floor() as i32 );
            let (frac_x, frac_y) = (x_pre.fract(), y_pre.fract());
            let v0 = self.dye_field[coord_to_idx((qx, qy))];
            let v1 = self.dye_field[coord_to_idx((qx + 1, qy))];
            let v2 = self.dye_field[coord_to_idx((qx, qy + 1))];
            let v3 = self.dye_field[coord_to_idx((qx + 1, qy + 1))];
            updated_dye_field[i] = blerp(v0, v1, v2, v3, frac_x, frac_y);
        }
        self.dye_field = updated_dye_field;
    }

    fn diffuse_velocities(&mut self, dt: f32) { // may be wrong
        let rate = dt * VISC_THINGY;// * (SIZE as f32 * SIZE as f32);
        for k in 0..GAUSS_SEIDEL_STEPS { // is 20 in paper
            let mut updated_vel_field = vec![(0.0,0.0); FLAT_SIZE];
            for (i, dye) in self.dye_field.iter().enumerate() {
                let (x, y) = idx_to_coord(i);
                updated_vel_field[i].0 = (
                    self.velocity_field[i].0 + rate * (
                        self.velocity_field[coord_to_idx((x as i32   , y as i32 +1))].0 +
                        self.velocity_field[coord_to_idx((x as i32   , y as i32 -1))].0 +
                        self.velocity_field[coord_to_idx((x as i32 +1, y as i32   ))].0 +
                        self.velocity_field[coord_to_idx((x as i32 -1, y as i32   ))].0
                    )) / (1.0 + 4.0 * rate);

                updated_vel_field[i].1 = (
                    self.velocity_field[i].1 + rate * (
                        self.velocity_field[coord_to_idx((x as i32   , y as i32 +1))].1 +
                        self.velocity_field[coord_to_idx((x as i32   , y as i32 -1))].1 +
                        self.velocity_field[coord_to_idx((x as i32 +1, y as i32   ))].1 +
                        self.velocity_field[coord_to_idx((x as i32 -1, y as i32   ))].1
                    )) / (1.0 + 4.0 * rate);
            }
            self.velocity_field = updated_vel_field;
        }
    }

    fn advect_velocities(&mut self, dt: f32) {
        let mut updated_vel_field = vec![(0.0,0.0); FLAT_SIZE];
        for (i, (vx, vy)) in self.velocity_field.iter().enumerate() {
            let (x, y) = idx_to_coord(i);
            let (x_pre, y_pre) = (
                (x as f32 - dt * vx).rem_euclid(SIZE as f32),
                (y as f32 - dt * vy).rem_euclid(SIZE as f32));
            let (qx, qy) = (x_pre.floor() as i32 , y_pre.floor() as i32 );
            let (frac_x, frac_y) = (x_pre.fract(), y_pre.fract());
            let (v0x, v0y) = self.velocity_field[coord_to_idx((qx, qy))];
            let (v1x, v1y) = self.velocity_field[coord_to_idx((qx + 1, qy))];
            let (v2x, v2y) = self.velocity_field[coord_to_idx((qx, qy + 1))];
            let (v3x, v3y) = self.velocity_field[coord_to_idx((qx + 1, qy + 1))];
            updated_vel_field[i] = (blerp(v0x, v1x, v2x, v3x, frac_x, frac_y),
                                    blerp(v0y, v1y, v2y, v3y, frac_x, frac_y));
        }
        self.velocity_field = updated_vel_field;
    }

    fn remove_divergence(&mut self) {
        let mut p = vec![0.0; FLAT_SIZE];
        let h = 1.0 / SIZE as f32;

        // find divergence
        let mut div = vec![0.0; FLAT_SIZE];
        for i in 0..FLAT_SIZE {
            let (x, y) = idx_to_coord(i);
            let (x, y) = (x as i32, y as i32);
            
            let ni = coord_to_idx((x, y + 1));
            let si = coord_to_idx((x, y - 1));
            let ei = coord_to_idx((x + 1, y));
            let wi = coord_to_idx((x - 1, y));

            let nv = self.velocity_field[ni];
            let sv = self.velocity_field[si];
            let ev = self.velocity_field[ei];
            let wv = self.velocity_field[wi];
  
            div[i] = -0.5 * h * ((nv.1 - sv.1) + (ev.1 - wv.1)); // might be rotated idk
        }

        // solve poisson
        for k in 0..GAUSS_SEIDEL_STEPS {
            for i in 0..FLAT_SIZE {
                let (x, y) = idx_to_coord(i);
                let (x, y) = (x as i32, y as i32);
                let n = p[coord_to_idx((x, y-1))];
                let s = p[coord_to_idx((x, y+1))];
                let e = p[coord_to_idx((x+1, y))];
                let w = p[coord_to_idx((x-1, y))];
                p[i] = (div[i] + w + e + n + s) / 4.0;
            }
        }

        for i in 0..FLAT_SIZE {
            let (x, y) = idx_to_coord(i);
            let (x, y) = (x as i32, y as i32);

            let n = p[coord_to_idx((x, y-1))];
            let s = p[coord_to_idx((x, y+1))];
            let e = p[coord_to_idx((x+1, y))];
            let w = p[coord_to_idx((x-1, y))];

            self.velocity_field[i].0 -= 0.5 * (e + w) / SIZE as f32;
            self.velocity_field[i].1 -= 0.5 * (n + s) / SIZE as f32
        }
    }

    fn update(&mut self) {
        let delta_time = self.current_time.elapsed();
        let dt = delta_time.as_secs_f32();
        
        self.diffuse_velocities(dt);    // WRONG

        // remove divergence before and after advection of velocities
        self.remove_divergence();       // WRONG
        self.advect_velocities(dt);     // WRONG
        self.remove_divergence();       // WRONG

        self.diffuse_dye(dt);
        self.advect_dye(dt);

        if delta_time < MIN_DT { sleep(MIN_DT - delta_time); }
        self.current_time = Instant::now();
    }

    fn events(&mut self) {
        for event in self.event_pump.poll_iter() {
            match event {
                Event::Quit {..} => { self.running = false; },
                _ => {}
            }
        }
        let ms = self.event_pump.mouse_state();
        let (x, y) = (ms.x() / RATIO as i32, ms.y() / RATIO as i32);
        if ms.left() {
            self.dye_field[coord_to_idx((x, y))] = 1.0;
        }
        if ms.right() {
            let rms = self.event_pump.relative_mouse_state();
            let (rx, ry) = (rms.x(), rms.y());
            self.velocity_field[coord_to_idx((x, y))].0 += rx as f32 * 1_000_000.0; // TODO: Why so small?
            self.velocity_field[coord_to_idx((x, y))].1 += ry as f32 * 1_000_000.0;
        }
    }

    fn render(&mut self) {
        self.canvas.set_draw_color(Color::BLACK);
        self.canvas.clear();
        for (i, (dye, (vx, vy))) in self.dye_field.iter().zip(self.velocity_field.iter()).enumerate() {
            if *dye > 0.0 {
                let (x, y) = field_to_canvas_corner(idx_to_coord(i));  
                let clr = (dye.powf(1./3.) * 255.0) as u8;
                self.canvas.set_draw_color(Color::RGB(0, clr, 0));
                self.canvas.fill_rect(Rect::from((x, y, RATIO, RATIO))).unwrap();
            }
            let (x, y) = idx_to_coord(i);
            let (x, y) = field_to_canvas_center(idx_to_coord(i));  
            let len = 1.0 + (vx * vx + vy * vy).sqrt(); // +1 to prevent divide by zero due to f32 imprecisions
            let (vx, vy) = ((6.0 * vx / len), (6.0 * vy / len));
            self.canvas.set_draw_color(Color::RGB(len as u8, 0, 255 - len as u8));
            self.canvas.draw_line((x, y), (x + vx as i32, y + vy as i32)).unwrap();
        }
        self.canvas.present();
    }
}

/* Idea?
 * Combine fluid solver (with two way coupling) with rigid and soft body solver
 */

/* Fluid solver
 * 1. Diffuse velocities
 * 2. Advect velocities
 * 3. Remove divergence
 * 
 * n. Diffuse dye
 * n+1. Advect die
 * Dye is for visualization purposes only, so doesnt really matter once boundary simulation is fulfilled.
*/