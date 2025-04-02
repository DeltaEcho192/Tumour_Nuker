use std::thread;
fn loop_test() {
    println!("Run Loop test");
    const grid_x: i32 = 2000;
    const grid_y: i32 = 4000;
    const grid_z: i32 = 100;
    const total_size: usize = (grid_x * grid_y * grid_z) as usize;

    let mut inbeam_count: i64 = 0;
    let mut operationcounter: i64 = 0;
    let mut handle_vec = Vec::new();
    let beam_amt: i32 = 12;
    for i in 1..beam_amt{
        handle_vec.push(thread::spawn( || {
            let mut dose: Vec<f32> = Vec::with_capacity(total_size);
            let beam_vec: Vector = Vector::new(0.8, 1.43, 6.3);

            let beam_ent_x: i32 = 2;
            let beam_ent_y: i32 = 2;
            let beam_ent_z: i32 = 0;

            let self_dot = beam_vec.dot(&beam_vec);

            let mut inbeam_count: i64 = 0;
            let mut operationcounter: i64 = 0;
            for x in 0..grid_x {
                for y in 0..grid_y {
                    for z in 0..grid_z {
                        let mut dose_val: f32 = 0.0;
                        let mut vector = Vector::new(x as f32,y as f32,z as f32);
                        vector.calculate_offset(beam_ent_x, beam_ent_y, beam_ent_z);
                        let dist = vector.dist_to_beam();
                        let dot_prod = vector.dot(&beam_vec);
                        let projection_point = vector.mult_vec(dot_prod/self_dot);
                        let project_dist = vector.dist_to_vector(&projection_point);
                        if project_dist <= 0.75 {
                            inbeam_count += 1; 
                            dose_val = project_dist;
                        }
                        //let index: usize = (x + grid_y * (y + grid_z * z)) as usize;
                        dose.push(dose_val);
                        operationcounter += 1;
                    }
                }
            }
        }));
    }

    for handle in handle_vec {
        handle.join().unwrap();
    }
    //}
//    println!("End Loop test");
//    println!("Operation Counter: {}", operationcounter);
//    let mut counter:i32 = 0;
//    for dose_val in &dose {
//        if *dose_val != 0.0 {
//            counter += 1;
//        }
//    }
//
//    println!("Counter {}", counter);
//    println!("Length {}", dose.len());
}


#[derive(Debug, Clone)]
struct Vector {
    x: f32,
    y: f32,
    z: f32,
}

impl Vector {
    fn new(x_in: f32, y_in: f32, z_in: f32) -> Vector {
        Vector {x: x_in, y: y_in, z: z_in}
    }
    fn dot(&self, v2: &Vector) -> f32 {
        self.x * v2.x + self.y * v2.y + self.z * v2.z
    }
    fn calculate_offset(&mut self, beam_x: i32, beam_y: i32, beam_z: i32) {
        self.x =  self.x - beam_x as f32;
        self.y =  self.y - beam_y as f32;
        self.z =  self.z - beam_z as f32;
    }
    fn dist_to_beam(&self) -> f32 {
        let dist:f32 = (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)) as f32;
        dist.sqrt()
    }

    fn dist_to_vector(&self, v2: &Vector) -> f32 {
        let dist:f32 = ((v2.x - self.x).powf(2.0) + (v2.y - self.y).powf(2.0) + (v2.z - self.z).powf(2.0)) as f32;
        dist.sqrt()
    }
    fn mult_vec(&mut self, val: f32) -> Vector {
        Vector {
        x : self.x * val,
        y : self.y * val,
        z : self.z * val,
    }
    }
}

fn main() {
    loop_test();
    println!("Hello, world!");
}
