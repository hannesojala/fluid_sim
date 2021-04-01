pub mod fluid;
use fluid::{*};
use std::{ thread::sleep, time::{Duration, Instant} };
use sdl2::{ event::Event, keyboard::Keycode, pixels::Color, rect::Rect };

const VISC: f32 = 1e-5;
const DIFF: f32 = 1e-5;

const SCALE: i32 = 3;
const SIZE: i32 = 256;
const MAX_FPS: u64 = 60;

pub struct Engine {
    pub running: bool,
    paused: bool,
    canvas: sdl2::render::WindowCanvas,
    event_pump: sdl2::EventPump,
    time: Instant,
    dt_s: f32,
    fluid: Fluid,
    draw_mode: i32
}

impl Engine {
    pub fn init() -> Engine {
        let sdl_context = sdl2::init().unwrap();
        let video_subsystem = sdl_context.video().unwrap();
        let win_size = (SIZE * SCALE) as u32;
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
            fluid: Fluid::new(VISC, DIFF, SIZE),
            draw_mode: 0
        }
    }

    // Updates the simulation state
    pub fn update(&mut self) {
        // Get time elapsed for timestep
        let delta_time = self.time.elapsed();
        println!("{}ms", delta_time.as_millis());
        self.time = Instant::now();
        if !self.paused {
            self.dt_s = delta_time.as_secs_f32(); // needs to be set for some ui funcs
            self.fluid.update(self.dt_s);
        }
        // If frame was faster than max framerate, sleep a lil
        static TARGET_DT : Duration = Duration::from_millis(1000 / MAX_FPS);
        if delta_time < TARGET_DT { sleep(TARGET_DT - delta_time); }
    }

    // Handles user input
    pub fn events(&mut self) {
        // Pump events in queue
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

        // Get both relative and global mouse coords
        let ms = self.event_pump.mouse_state();
        let (mx, my) = (ms.x() / SCALE , ms.y() / SCALE);
        let rms = self.event_pump.relative_mouse_state();
        let (rx, ry) = (rms.x(), rms.y());

        // Get box for draw tool radius
        let radius: i32 = SIZE / 32;
        let (y0,y1) = ((my-radius).max(1), (my+radius).min(SIZE-2));
        let (x0,x1) = ((mx-radius).max(1), (mx+radius).min(SIZE-2));
        
        // For all in draw tool box
        for y in y0..y1 {
            for x in x0..x1 {
                // If dist <= radius
                if (((mx-x)*(mx-x)+(my-y)*(my-y))as f32).sqrt() <= radius as f32 {
                    // Add dye or velocity
                    if ms.left() {
                        self.fluid.add_vel(x, y, self.dt_s * (rx * SIZE) as f32, 
                                                 self.dt_s * (ry * SIZE) as f32);
                    }
                    if ms.right() {
                        self.fluid.add_dye(x, y, self.dt_s * 32.0);
                    }
                }
            }
        }
    }

    // Renders either vector field or dye field inefficiently
    // TODO: This is a terrible way to do this.
    pub fn render(&mut self) {
        self.canvas.set_draw_color(Color::BLACK);
        self.canvas.clear();
        for y in 0..SIZE {
            for x in 0..SIZE {
                // Draw dye field density
                let dye_amt = self.fluid.dye_at(x,y);
                if  self.draw_mode == 1 && dye_amt > 0.0 {
                    let color = Color::RGB(0, 0, (dye_amt*255.) as u8); // sqrt() because of percieved brightness (gamma)
                    let rect = Rect::new(SCALE * x, SCALE * y, SCALE as u32, SCALE as u32);
                    self.canvas.set_draw_color(color);
                    self.canvas.fill_rect(rect).unwrap();
                }
                // Draw vector field component intensity
                let (vx, vy) = self.fluid.vel_at(x,y);
                if self.draw_mode == 0 && (vx != 0. || vy != 0.) {
                    let color = Color::RGB( (255.0 * vx / SIZE as f32).abs() as u8, // abs() for either direction
                                            (255.0 * vy / SIZE as f32).abs() as u8, 0 );
                    let rect = Rect::new(SCALE * x, SCALE * y, SCALE as u32, SCALE as u32);
                    self.canvas.set_draw_color(color);
                    self.canvas.fill_rect(rect).unwrap();
                }
            }
        }
        self.canvas.present();
    }
}