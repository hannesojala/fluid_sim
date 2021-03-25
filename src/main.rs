/* TODO
 * - Set edges to be fluid boundaries, see if that improves simulation -> IT SHOULD I THINK
 * - Image Dye
 * - generalize gauss-seidel code
 * - generalize diffuse code
 * - generalize advect code
 * - check for mistakes
*/

use std::thread::sleep;
use std::time::{Duration, Instant};
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Rect;

const N: i32 = 256;
const SIZE: usize = (N*N) as usize;
const RATIO: i32 = 2;
const MAX_FPS: u64 = 1000;
const MIN_DT : Duration = Duration::from_millis(1000 / MAX_FPS);

const VISC: f32 = 5e-2;
const DIFF: f32 = 5e-9;
const APPROX: u32 = 20;
const STIR: f32 = 1e-4;

fn idx(x: i32, y: i32) -> usize { (N * y.rem_euclid(N) + x.rem_euclid(N)) as usize }

fn coord(i: usize) -> (i32, i32) { (i as i32 % N, i as i32 / N) }

fn lerp(v1: f32, v2: f32, k: f32) -> f32 { v1 + k * (v2 - v1) }

fn blerp(v1: f32, v2: f32, v3: f32, v4: f32, k1: f32, k2: f32) -> f32 { 
    lerp(lerp(v1, v2, k1), lerp(v3, v4, k1), k2)
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
    time: Instant,
    dye: Vec<f32>,
    vel: Vec<(f32, f32)>,
    draw_mode: i32
}

impl Engine {
    fn init() -> Engine {
        let sdl_context = sdl2::init().unwrap();
        let video_subsystem = sdl_context.video().unwrap();
        let win_size = (N * RATIO) as u32;
        let window = video_subsystem.window("fluid", win_size, win_size)
            .position_centered()
            .build().unwrap();
        
        Engine {
            canvas: window.into_canvas().build().unwrap(),
            event_pump: sdl_context.event_pump().unwrap(),
            running: true,
            time: Instant::now(),
            dye: vec![0.0; SIZE],
            vel: vec![(0.0,0.0); SIZE],
            draw_mode: 0
        }
    }

    fn diffuse_dye(&mut self, dt: f32) {
        let rate = dt * DIFF * (N*N) as f32;
        for _k in 0..APPROX {
            let mut updated_dye_field = vec![0.0; SIZE]; // TODO DONT REALLOC
            for y in 0..N {
                for x in 0..N {
                    updated_dye_field[idx(x, y)] = (self.dye[idx(x, y)] + rate * (
                        self.dye[idx(x  , y+1)] +
                        self.dye[idx(x  , y-1)] +
                        self.dye[idx(x+1, y  )] +
                        self.dye[idx(x-1, y  )]
                    )) / (1.0 + 4.0 * rate);
                }
            }
            self.dye = updated_dye_field;
        }
    }

    fn advect_dye(&mut self, dt: f32) {
        let mut updated_dye_field = vec![0.0; SIZE];
        for y in 0..N {
            for x in 0..N {
                let (vx, vy) = self.vel[idx(x,y)];
                let (x0, y0) = (
                    (x as f32 - dt * vx).rem_euclid(N as f32),
                    (y as f32 - dt * vy).rem_euclid(N as f32));
                let (qx, qy) = (x0.floor() as i32, y0.floor() as i32);
                let v0 = self.dye[idx(qx  , qy  )];
                let v1 = self.dye[idx(qx+1, qy  )];
                let v2 = self.dye[idx(qx  , qy+1)];
                let v3 = self.dye[idx(qx+1, qy+1)];
                updated_dye_field[idx(x ,y)] = blerp(v0, v1, v2, v3, x0.fract(), y0.fract());
            }
        }
        self.dye = updated_dye_field;
    }

    fn diffuse_velocities(&mut self, dt: f32) {
        let rate = dt * VISC * (N*N) as f32;
        for _k in 0..APPROX {
            let mut updated_vel_field = vec![(0.0,0.0); SIZE];
            for y in 0..N {
                for x in 0..N {
                    updated_vel_field[idx(x,y)].0 = (
                        self.vel[idx(x,y)].0 + rate * (
                            self.vel[idx(x  , y+1)].0 +
                            self.vel[idx(x  , y-1)].0 +
                            self.vel[idx(x+1, y  )].0 +
                            self.vel[idx(x-1, y  )].0
                        )) / (1.0 + 4.0 * rate);

                    updated_vel_field[idx(x,y)].1 = (
                        self.vel[idx(x,y)].1 + rate * (
                            self.vel[idx(x  , y+1)].1 +
                            self.vel[idx(x  , y-1)].1 +
                            self.vel[idx(x+1, y  )].1 +
                            self.vel[idx(x-1, y  )].1
                        )) / (1.0 + 4.0 * rate);
                    }
            }
            self.vel = updated_vel_field;
        }
    }

    fn advect_velocities(&mut self, dt: f32) {
        let mut updated_vel_field = vec![(0.0,0.0); SIZE];
        for y in 0..N {
            for x in 0..N {
                let (vx, vy) = self.vel[idx(x,y)];
                let (x_pre, y_pre) = (
                    (x as f32 - (dt*N as f32) * vx).rem_euclid(N as f32), // mul dt by N according to stam, maybe not needed if boxed in fluid?
                    (y as f32 - (dt*N as f32) * vy).rem_euclid(N as f32)); // instead of remainder he corrects (to 0.5 and n-0.5)
                let (qx, qy) = (x_pre.floor() as i32, y_pre.floor() as i32);
                let (frac_x, frac_y) = (x_pre.fract(), y_pre.fract());
                let (v0x, v0y) = self.vel[idx(qx  , qy  )];
                let (v1x, v1y) = self.vel[idx(qx+1, qy  )];
                let (v2x, v2y) = self.vel[idx(qx  , qy+1)];
                let (v3x, v3y) = self.vel[idx(qx+1, qy+1)];
                updated_vel_field[idx(x,y)] = (blerp(v0x, v1x, v2x, v3x, frac_x, frac_y), // can be condensed i think (see stam)
                                        blerp(v0y, v1y, v2y, v3y, frac_x, frac_y));
            }
        }
        self.vel = updated_vel_field;
    }

    fn remove_divergence(&mut self) {
        let mut p = vec![0.0; SIZE];
        let mut div = vec![0.0; SIZE];
        for y in 0..N {
            for x in 0..N {
                let nv = self.vel[idx(x, y + 1)];
                let sv = self.vel[idx(x, y - 1)];
                let ev = self.vel[idx(x + 1, y)];
                let wv = self.vel[idx(x - 1, y)];
                div[idx(x,y)] = -0.5 * (sv.1 - nv.1 + ev.0 - wv.0) / N as f32;
            }
        }

        for _k in 0..APPROX {
            for y in 0..N {
                for x in 0..N {
                    p[idx(x,y)] = (div[idx(x,y)] + p[idx(x, y-1)] + p[idx(x, y+1)]
                                                 + p[idx(x+1, y)] + p[idx(x-1, y)]) / 4.0;
                }
            }
        }

        for y in 0..N {
            for x in 0..N {
                self.vel[idx(x,y)].0 -= 0.5 * (p[idx(x+1, y)] + p[idx(x-1, y)]) * N as f32; // 4dummies says * N (same, same)
                self.vel[idx(x,y)].1 -= 0.5 * (p[idx(x, y-1)] + p[idx(x, y+1)]) * N as f32; // stam says / h, where h = 1/N, but now it explodes!
            }
        }
    }

    fn update(&mut self) {
        let delta_time = self.time.elapsed();
        self.time = Instant::now();
        println!("{}ms", delta_time.as_millis());
        let dt = delta_time.as_secs_f32();
        
        self.dye[idx(N/2, N/2)] = 1.0;

        self.diffuse_velocities(dt);    
        self.remove_divergence();       
        self.advect_velocities(dt);     
        self.remove_divergence();       

        self.diffuse_dye(dt);
        self.advect_dye(dt);

        if delta_time < MIN_DT { sleep(MIN_DT - delta_time); }
    }

    fn events(&mut self) {
        for event in self.event_pump.poll_iter() {
            match event {
                Event::Quit {..} => { self.running = false; },
                Event::KeyDown { keycode: Some(Keycode::M), repeat: false, .. } => {
                    self.draw_mode = (self.draw_mode + 1) % 2;
                },
                _ => {}
            }
        }
        let ms = self.event_pump.mouse_state();
        let (mx, my) = (ms.x() / RATIO , ms.y() / RATIO);
        if ms.left() {
            let rms = self.event_pump.relative_mouse_state();
            let (rx, ry) = (rms.x(), rms.y());
            self.vel[idx(mx, my)].0 += rx as f32 * STIR; // TODO: Why so small?
            self.vel[idx(mx, my)].1 += ry as f32 * STIR
        }
    }

    fn render(&mut self) {
        self.canvas.set_draw_color(Color::BLACK);
        self.canvas.clear();
        for (i, (dye, (vx, vy))) in self.dye.iter().zip(self.vel.iter()).enumerate() {
            let (x, y) = coord(i);  
            if  self.draw_mode == 1 && *dye > 0.0 {
                self.canvas.set_draw_color(Color::RGB(0, (*dye * 255.0) as u8, 0));
                self.canvas.fill_rect(Rect::from((RATIO * x, RATIO * y, RATIO as u32, RATIO as u32))).unwrap();
            }
            if self.draw_mode == 0 && (*vx != 0. || *vy != 0.) {
                self.canvas.set_draw_color(Color::RGB(vx.abs().min(255.) as u8, vy.abs().min(255.) as u8, 0));
                self.canvas.fill_rect(Rect::from((RATIO * x, RATIO * y, RATIO as u32, RATIO as u32))).unwrap();
            }
        }
        self.canvas.present();
    }
}