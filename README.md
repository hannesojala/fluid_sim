# fluid simulator
A fluid simulator designed to learn Rust.

Switch between viewing modes (fluid velocity, default=dye color) by pressing "M".
Switch between dye colors (default=Red, Green, Blue and Eraser Dye) by pressing "C".
Pause with 'Pause'.

If an image called 'image.jpg' is in the executables directory it will be loaded as the background. This is implemented very lazily, I know.

Swirl the liquid with mouse-left, and add dye of the selected color with mouse-right.

![alt text](demo.png?raw=true)
![alt text](demo2.png?raw=true)

To compile: Install Rust, the easiest way is with [rustup](https://rustup.rs/), and simply build and run with the command "cargo run --release". If your OS does not have the required SDL2 libraries, follow the SDL2 install instructions here: https://github.com/Rust-SDL2/rust-sdl2. On windows you might need to add a SDL.dll file. Additionally, using the "bundled" feature described in the previous link will probably alleviate all of this trouble.
