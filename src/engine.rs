#[macro_use]
pub mod util;
use util::{*};
use std::{ thread::sleep, time::{Duration, Instant} };
use sdl2::{ event::Event, keyboard::Keycode, pixels::Color, rect::Rect };

pub struct Engine {
    pub running: bool,
    paused: bool,
    canvas: sdl2::render::WindowCanvas,
    event_pump: sdl2::EventPump,
    time: Instant,
    dt_s: f32,
    dye : Vec<f32>,
    vel : Vec<(f32, f32)>,
    draw_mode: i32
}

impl Engine {
    pub fn init() -> Engine {
        let sdl_context = sdl2::init().unwrap();
        let video_subsystem = sdl_context.video().unwrap();
        let win_size = (N * RATIO) as u32;
        let window = video_subsystem.window("fluid", win_size, win_size)
            .position_centered()
            .build().unwrap();
        
        Engine {
            running: true,
            paused: false,
            canvas: window.into_canvas().build().unwrap(),
            event_pump: sdl_context.event_pump().unwrap(),
            time: Instant::now(),
            dt_s: 0.0,
            dye:  vec![0.0; SIZE],
            vel:  vec![(0.0,0.0); SIZE],
            draw_mode: 0
        }
    }

    pub fn update(&mut self) {
        let delta_time = self.time.elapsed();
        self.time = Instant::now();
        if !self.paused {
            self.dt_s = delta_time.as_secs_f32();

            // let mass_before: f32 = self.dye.iter().sum();

            /* DIFFUSE VELOCITY FIELD*/
            self.vel = solve_vfield(&self.vel, self.dt_s * VISC * ((N - 2)*(N - 2)) as f32);
            // remove div
            self.vel =  enforce_div_eq_0(&self.vel);

            /* ADVECT VELOCITY FIELD*/
            self.vel = advect_vfield(&self.vel, self.dt_s);
            // remove div
            self.vel =  enforce_div_eq_0(&self.vel);

            /* DIFFUSE DYE FIELD*/
            self.dye = solve_sfield(&self.dye, self.dt_s * DIFF * ((N - 2)*(N - 2)) as f32);

            /* ADVECT DYE FIELD */
            self.dye = advect_sfield(&self.dye, &self.vel, self.dt_s);


            // let mass_after: f32 = self.dye.iter().sum();
            // fails sometimes due to advection
            // works with diffusion
            // todo: enforce mass conservation of dye during advection
            // assert!((mass_before - mass_after).abs() <= 0.01 * mass_before); // Make sure dye mass is conserved

        }
        static TARGET_DT : Duration = Duration::from_millis(1000 / MAX_FPS);
        if delta_time < TARGET_DT { sleep(TARGET_DT - delta_time); }
    }

    pub fn events(&mut self) {
        for event in self.event_pump.poll_iter() {
            match event {
                Event::Quit {..} => { self.running = false; },
                Event::KeyDown { keycode: Some(Keycode::M), repeat: false, .. } => {
                    self.draw_mode = (self.draw_mode + 1) % 2;
                },
                Event::KeyDown { keycode: Some(Keycode::Pause), repeat: false, .. } => {
                    self.paused = !self.paused;
                },
                _ => {}
            }
        }

        let ms = self.event_pump.mouse_state();
        let (mx, my) = (ms.x() / RATIO , ms.y() / RATIO);
        let rms = self.event_pump.relative_mouse_state();
        let (rx, ry) = (rms.x(), rms.y());

        let radius: i32 = 7;
        let (y0,y1) = ((my-radius).max(1), (my+radius).min(N-2));
        let (x0,x1) = ((mx-radius).max(1), (mx+radius).min(N-2));
        
        for y in y0..y1 { // restrict to bounds
            for x in x0..x1 {
                if (((mx-x)*(mx-x)+(my-y)*(my-y))as f32).sqrt() <= radius as f32 { // radius
                    if ms.left() {
                        self.vel[i!(x, y)].0 += self.dt_s * (rx * N / RATIO) as f32; // Divide by ratio to get
                        self.vel[i!(x, y)].1 += self.dt_s * (ry * N / RATIO) as f32; // cells traversed by mouse
                    }
                    if ms.right() {
                        self.dye[i!(x, y)] += self.dt_s * 16.0; // should make dt based to account for frame rate and pause
                    }
                }
            }
        }
    }

    pub fn render(&mut self) {
        self.canvas.set_draw_color(Color::BLACK);
        self.canvas.clear();
        for y in 0..N {
            for x in 0..N {
                let dye_amt = self.dye[i!(x,y)];
                if  self.draw_mode == 1 && dye_amt > 0.0 {
                    let color = Color::RGB(0, (dye_amt.sqrt()*255.) as u8, 0); // dont need dye_amt.sqrt() if add more dye
                    let rect = Rect::new(RATIO * x, RATIO * y, RATIO as u32, RATIO as u32);
                    self.canvas.set_draw_color(color);
                    self.canvas.fill_rect(rect).unwrap();
                }
                let (vx, vy) = self.vel[i!(x,y)];
                if self.draw_mode == 0 && (vx != 0. || vy != 0.) {
                    let color = Color::RGB( (255.0 * vx / N as f32).abs() as u8,         // * 8 / N to get ratio of vel
                                            (255.0 * vy / N as f32).abs() as u8, 0 );    // cmp'd to vel from 8 cell mouse drag
                    let rect = Rect::new(RATIO * x, RATIO * y, RATIO as u32, RATIO as u32);
                    self.canvas.set_draw_color(color);
                    self.canvas.fill_rect(rect).unwrap();
                }
            }
        }
        self.canvas.present();
    }
}