mod engine;
fn main() {
    let mut engine = engine::Engine::init();
    while engine.running {
        engine.events();
        engine.update();
        engine.render();
    }
}