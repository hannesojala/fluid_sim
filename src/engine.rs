pub mod fluid;
use fluid::{*};
use image::{GenericImageView, imageops::FilterType::Nearest};
use std::{thread::sleep, time::{Duration, Instant}};
use sdl2::{event::Event, keyboard::Keycode, pixels::Color, rect::Rect};

const VISC: f32 = 1e-15;
const DIFF: f32 = 1e-15;
const VORT: f32 = 5.;

const SCALE: i32 = 3;
const SIZE: i32 = 256;
const MAX_FPS: u64 = 144;

// Quality related, 5 for low, 10-15 for medium, 20 for high
const GAUSS_SEIDEL_ITERATIONS: u32 = 15;

const COLORS: [(i16, i16, i16); 4] = [
    (255, 0, 0),
    (0, 255, 0),
    (0, 0, 255),
    (-255, -255, -255),
];

pub struct Engine {
    pub running: bool,
    paused: bool,
    canvas: sdl2::render::WindowCanvas,
    event_pump: sdl2::EventPump,
    fluid: Fluid,
    draw_mode: i32,
    time_frame_start: Instant,
    delta_time: Duration,
    total_time: Duration,
    draw_color_index: usize,
}

impl Engine {
    pub fn init() -> Engine {
        let sdl_context = sdl2::init().unwrap();
        let video_subsystem = sdl_context.video().expect("Put an image of size SIZE in this directory called image.jpg");
        let win_size = (SIZE * SCALE) as u32;
        let window = video_subsystem.window("fluid", win_size, win_size)
            .position_centered()
            .build().unwrap();
        let canvas = window.into_canvas().build().unwrap();
        let im = image::open(&std::path::Path::new("./image.jpg"));
        let mut fluid = Fluid::new(VISC, DIFF, VORT, GAUSS_SEIDEL_ITERATIONS, SIZE);
        if let Ok(i) = im {
            let image = i.resize_to_fill(SIZE as u32, SIZE as u32, Nearest);
            for y in 0..SIZE {
                for x in 0..SIZE {
                    let pixel = image.get_pixel(x as u32, y as u32);
                    fluid.set_dye(x, y, (pixel.0[0] as f32, pixel.0[1] as f32, pixel.0[2] as f32 ))
                }
            }
        }   
        
        Engine {
            running: true,
            paused: false,
            canvas,
            event_pump: sdl_context.event_pump().expect("Could not get event pump!"),
            fluid,
            draw_mode:          1,
            time_frame_start:   Instant::now(),
            delta_time:         Duration::from_millis(0),
            total_time:         Duration::from_millis(0),
            draw_color_index:   0
        }
    }

    // Updates the simulation state WRONG TIMESTEP CODE ITS BROKEN
    pub fn update(&mut self) {
        // Get time elapsed for timestep
        self.time_frame_start = Instant::now();
        self.total_time += self.delta_time;
        println!("{}ms", self.delta_time.as_millis());
        if !self.paused {
            self.fluid.update(self.delta_time.as_secs_f32());
        }
        static TARGET_DT: Duration = Duration::from_millis(1000 / MAX_FPS);
        let frame_time = Instant::now() - self.time_frame_start;
        if frame_time < TARGET_DT {
            sleep(TARGET_DT - frame_time);
        }
        self.delta_time = Instant::now() - self.time_frame_start;
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
                Event::KeyDown { keycode: Some(Keycode::C), repeat: false, .. } => {
                    self.draw_color_index = (self.draw_color_index + 1) % COLORS.len();
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
        
        let dt_s = self.delta_time.as_secs_f32();
        
        // For all in draw tool box
        for y in y0..y1 {
            for x in x0..x1 {
                // If dist <= radius
                if (((mx-x)*(mx-x)+(my-y)*(my-y))as f32).sqrt() < radius as f32 {
                    // Add dye or velocity
                    if ms.left() {
                        self.fluid.add_vel(x, y, dt_s * ((mx-x).abs() * rx * SIZE) as f32, 
                                                 dt_s * ((my-y).abs() * ry * SIZE) as f32);
                    }
                    if ms.right() {
                        let clr = COLORS[self.draw_color_index];
                        let dye = (
                            clr.0 as f32,
                            clr.1 as f32,
                            clr.2 as f32,
                        );
                        self.fluid.set_dye(x, y, dye);
                    }
                    if ms.middle() {
                        
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
                let dye = self.fluid.dye_at(x,y);
                if  self.draw_mode == 1 {
                    let color = Color::RGB(dye.0 as u8, dye.1 as u8, dye.2 as u8);
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